use std::collections::HashMap;
use crate::errors::GenerateReadsErrors;
use crate::common::structs::{
    variants::{
        Variant, 
        VariantType, 
        Genotypes::{Homozygous, Heterozygous}
    },
    distributions::DiscreteDistribution,
};
use simple_rng::NeatRng;



#[cfg(test)]
mod tests {
    use super::*;

   
}