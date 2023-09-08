use rand::seq::SliceRandom;
use rand::thread_rng;

fn find_random_non_n(positions: Vec<u64>) -> u64 {
    /*
    This function is intended to find a non-n position from a vector
    An improvement is would be to use the gc-bias as a weighting vector, path them both in
    and then use weighted sample.
     */
    let mut rng = thread_rng();
    let selected = positions.choose(&mut rng);
    *selected.unwrap()
}

fn map_non_n_regions(sequence: String) -> Vec<bool> {
    let check_vec: Vec<bool> = vec![true; sequence.len()];
    for (index, char) in sequence.chars().enumerate() {
        if char == "N" {
            check_vec[index] = false;
        }
    }
    
    check_vec
}

fn main() {
    /*
    Ultimately, this vector will be a list of positions to sample from
     */
    let my_positions: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7];
    let a_position: u64 = find_random_non_n(my_positions);
    println!("Insert at {:?}", a_position);
}
