use btree_range_map;
use btree_range_map::Measure;

pub fn main() {
    btree_range_map::RangeMap::<i32, &str>::new();
    println!("{}", &<i32 as Measure<i32>>::len(&5));
    println!("{}", &<i32 as Measure<i32>>::distance(&3, &7));
}
