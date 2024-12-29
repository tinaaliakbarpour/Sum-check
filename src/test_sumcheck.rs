
use lambdaworks_math::{field::{element::FieldElement,
    fields::u64_prime_field::U64PrimeField}};


use crate::sumcheck;
use crate::polynomial::{Poly, Term};

const ORDER: u64 = 23;


type F = U64PrimeField<ORDER>;

#[test]
fn test_sumcheck() {
    // Constructing the two-variate multilinear polynomial x_0 + 2 * x_1 + 3 * x_0 * x_1
    let terms = vec![
        (FieldElement::<F>::one(), vec![Term{ index: 0, power:1}]), 
        (FieldElement::new(2), vec![Term{index:1, power:1}]),
        (FieldElement::new(3), vec![Term{index:0, power:1}, Term{index:1, power:1}]),
    ];
    let g = Poly::new(2,terms);

	sumcheck::verify(g, FieldElement::from(9));

}