use rand::seq::SliceRandom;
use rand::thread_rng;

fn find_random_non_n(positions: Vec<u64>) -> u64 {
    let mut rng = thread_rng();
    let selected = positions.choose(&mut rng);
    *selected.unwrap()
}

fn main() {
    let my_positions: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7];
    let a_position: u64 = find_random_non_n(my_positions);
    println!("Insert at {:?}", a_position);
}
