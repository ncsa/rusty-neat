//! Block-wise gzip writer backed by libdeflate.
//!
//! libdeflate compresses substantially faster than streaming zlib/zlib-ng for
//! the same gzip output, but its API is one-shot (compress a whole buffer),
//! not streaming. `BlockGzWriter` bridges that: it buffers incoming bytes and,
//! whenever it has accumulated ~`BLOCK_SIZE`, compresses that block into an
//! independent gzip *member* and writes it out. Concatenated gzip members form
//! a standard multi-member gzip stream — `gunzip`, `zcat`, `MultiGzDecoder`,
//! and every common bioinformatics tool read them transparently.
//!
//! Memory is bounded to one block regardless of total output size, preserving
//! rneat's streaming, low-memory behaviour while cutting compression CPU.

use libdeflater::{CompressionLvl, Compressor};
use log::error;
use std::io::{self, Write};

/// Uncompressed bytes buffered before a block is emitted. 256 KiB balances
/// compression ratio (bigger = slightly better) against memory and per-member
/// header overhead (smaller = more gzip headers).
const BLOCK_SIZE: usize = 256 * 1024;

/// A `Write` adapter that emits multi-member gzip via libdeflate.
///
/// Call [`BlockGzWriter::finish`] to flush the final partial block and observe
/// any error. As a safety net `Drop` also flushes (matching `flate2::GzEncoder`,
/// whose `Drop` likewise finalizes and swallows errors), so existing code that
/// relies on drop-finalization keeps working.
pub struct BlockGzWriter<W: Write> {
    inner: W,
    buf: Vec<u8>,
    compressor: Compressor,
    finished: bool,
    /// Whether any gzip member has been emitted. If none has by finalize time,
    /// we emit an empty member so the output is still a valid (empty) gzip
    /// stream — matching `GzEncoder`, and avoiding 0-byte files that a decoder
    /// would reject with "unexpected end of file" (e.g. a fully deleted contig).
    wrote_any: bool,
}

impl<W: Write> BlockGzWriter<W> {
    pub fn new(inner: W) -> Self {
        // CompressionLvl::default() is libdeflate level 6 — the same ratio
        // target as flate2's Compression::default(), so output size is
        // comparable and only the speed changes.
        Self {
            inner,
            buf: Vec::with_capacity(BLOCK_SIZE),
            compressor: Compressor::new(CompressionLvl::default()),
            finished: false,
            wrote_any: false,
        }
    }

    /// Compress `block` (which may be empty) as one gzip member and write it.
    fn emit_member(&mut self, block: &[u8]) -> io::Result<()> {
        let bound = self.compressor.gzip_compress_bound(block.len());
        let mut out = vec![0u8; bound];
        let n = self
            .compressor
            .gzip_compress(block, &mut out)
            .map_err(|e| io::Error::other(format!("libdeflate gzip compression failed: {e:?}")))?;
        self.inner.write_all(&out[..n])?;
        self.wrote_any = true;
        Ok(())
    }

    /// Emit the buffered remainder, and — if nothing was ever written — an empty
    /// member, so the stream is always valid gzip.
    fn flush_buffered(&mut self) -> io::Result<()> {
        if !self.buf.is_empty() {
            let block = std::mem::take(&mut self.buf);
            self.emit_member(&block)?;
        }
        if !self.wrote_any {
            self.emit_member(&[])?;
        }
        Ok(())
    }

    /// Flush the final partial block and the inner writer. After this, `Drop`
    /// is a no-op.
    pub fn finish(mut self) -> io::Result<()> {
        self.flush_buffered()?;
        self.finished = true;
        self.inner.flush()
        // self drops here; Drop sees `finished` and does nothing.
    }
}

impl<W: Write> Write for BlockGzWriter<W> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.buf.extend_from_slice(data);
        // Emit full blocks, keeping any remainder buffered.
        while self.buf.len() >= BLOCK_SIZE {
            let rest = self.buf.split_off(BLOCK_SIZE);
            let block = std::mem::replace(&mut self.buf, rest);
            self.emit_member(&block)?;
        }
        Ok(data.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        // Do not force a block boundary on every flush() — that would create
        // many tiny gzip members. Blocks are emitted at BLOCK_SIZE or finish().
        self.inner.flush()
    }
}

impl<W: Write> Drop for BlockGzWriter<W> {
    fn drop(&mut self) {
        if self.finished {
            return;
        }
        if let Err(e) = self.flush_buffered() {
            error!("BlockGzWriter: failed to flush final gzip block on drop: {e}");
        }
    }
}
