pub extern crate simplelog;
pub extern crate log;

pub use log::*;

pub use simplelog::{SimpleLogger, WriteLogger, LevelFilter, CombinedLogger, TerminalMode, Config, ColorChoice};

pub use std::fs::File;

pub fn init() {
    CombinedLogger::init(
        vec![
            #[cfg(feature = "termcolor")]
            TermLogger::new(
                LevelFilter::Warn,
                Config::default(),
                TerminalMode::Mixed,
                ColorChoice::Auto,
            ),
            #[cfg(not(feature="termolor"))]
            SimpleLogger::new(LevelFilter::Info, Config::default()),
            WriteLogger::new(
                LevelFilter::Debug,
                Config::default(),
                File::create("../this_log.log").unwrap()
            ),
        ]
    ).unwrap();
}
