pub struct ReadRecord {
    pub name: String,
    pub sequence: String,
    pub quality_scores: Vec<usize>,
    pub cigar_ops: Vec<char>,
    pub is_paired: bool,
    pub is_reverse: bool,
    pub contig: String,
    pub position: usize,
    pub mate_contig: String,
    pub mate_position: usize,
    pub template_length: i32,
}
