use core::num;
use std::vec;

use lambdaworks_math::{
    field::{element::FieldElement, traits::IsField},
};

use lambdaworks_math::polynomial::dense_multilinear_poly::DenseMultilinearPolynomial;
use lambdaworks_math::polynomial::Polynomial as UniPoly;
use rand::Rng;
use serde_json::value::Index;

use crate::polynomial::{Poly, Term};


// Converts i into an index in {0,1}^v
pub fn n_to_vec<F: IsField>(i: usize, n: usize) -> Vec<FieldElement<F>> {
	format!("{:0>width$}", format!("{:b}", i), width = n)
		.chars()
		.map(|x| if x == '1' { 1.into() } else { 0.into() })
		.collect()
}

// Simulates memory of a single prover instance
// #[derive(Debug, Clone)]
pub struct Prover <F: IsField>
where F::BaseType: Send + Sync{
	pub g: Poly<F>,
	pub r_vec: Vec<FieldElement<F>>,
}

impl <F: IsField > Prover<F>
where  
    F::BaseType: Send + Sync{
	pub fn new(g: Poly<F>) -> Self {
		Prover {
			g: g,
			r_vec: vec![],
		}
	}


	// Given polynomial g, fix Xj, evaluate over xj+1
	pub fn gen_uni_polynomial(&mut self, r: Option<FieldElement<F>>) -> UniPoly<FieldElement<F>> {
		if r.is_some() {
			self.r_vec.push(r.unwrap());
		}
		let v = self.g.num_vars() - self.r_vec.len();
        eprintln!("{:?}",v);
		(0..(2u32.pow(v as u32 - 1))).fold(
			UniPoly::zero(),
			|sum, n| sum + self.evaluate_gj(n_to_vec(n as usize, v)),
		)
	}
	// Evaluates gj over a vector permutation of points, folding all evaluated terms together into one univariate polynomial
	pub fn evaluate_gj(&self, points: Vec<FieldElement<F>>) -> UniPoly<FieldElement<F>> {
		self.g.terms.iter().fold(
			UniPoly::zero(),
			|sum, (coeff, term)| {
				let mut curr = self.evaluate_term(&term, &points);
                curr = curr * coeff;
				curr + sum
			},
		)
	}

	// Evaluates a term with a fixed univar, returning (new coefficent, fixed term)
	pub fn evaluate_term(
		&self,
		terms: &Vec<Term>,
		point: &Vec<FieldElement<F>>,
	) -> UniPoly<FieldElement<F>> {
		// let mut fixed_term: Option<SparseTerm> = None;
        let mut uni_term_power = 0 as usize; 
		let coeff: FieldElement<F> =
            terms.iter().fold(FieldElement::one(), |product, term| match term.index {
				j if j == self.r_vec.len() => {
					uni_term_power = term.power;
					product
				}
				j if j < self.r_vec.len() => self.r_vec[j].pow(term.power) * product,
				_ => point[term.index - self.r_vec.len()].pow(term.power) * product,
			});
		// (coeff, fixed_term)
        UniPoly::new_monomial(coeff, uni_term_power)
	}

}

// Verifier procedures
pub fn get_r<F: IsField>() -> Option<FieldElement<F>> {
	let mut rng = rand::thread_rng();
	let r: FieldElement<F> = FieldElement::from(rng.gen::<u64>());
    // let r = FieldElement::one();
	Some(r)
}

// Verify prover's claim c_1
// Presented pedantically:
pub fn verify<F: IsField>(g: Poly<F>, c_1: FieldElement<F>) -> bool
where F::BaseType: Send + Sync {
	// 1st round

	let mut p = Prover::new(g);
    let num_var = p.g.num_vars();
    // let mut gi  = p.g.clone();
    let mut gi = p.gen_uni_polynomial(None);
    eprintln!("{:?}",gi);
	let mut expected_c = gi.evaluate(&FieldElement::<F>::zero()) + gi.evaluate(&FieldElement::one()) ;
	assert_eq!(c_1, expected_c);
	

	// middle rounds
	for _ in 1..num_var {
		let r = get_r().unwrap();
		expected_c = gi.evaluate(&r);
        gi = p.gen_uni_polynomial(Some(r.clone()));
		let new_c = gi.evaluate(&FieldElement::<F>::one()) + gi.evaluate(&FieldElement::zero());
		assert_eq!(expected_c, new_c);
	}
	// final round
	let r = get_r().unwrap();

	expected_c = gi.evaluate(&r);
	p.r_vec.push(r);
	let new_c = p.g.evaluate(&p.r_vec);
	assert_eq!(expected_c, new_c);
	true
}

