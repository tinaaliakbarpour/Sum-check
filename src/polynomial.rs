use lambdaworks_math::{field::{element::FieldElement, traits::IsField}};

pub struct Term {
    pub index : usize,
    pub power:usize,
}
pub struct Poly<F: IsField> {
    pub max_degree: usize,
    pub terms: Vec<(FieldElement<F>,Vec<Term>)>,
    pub num_vars:usize,
}


impl <F: IsField>Poly<F> {

    pub fn new(num_vars: usize,terms: Vec<(FieldElement<F>, Vec<Term>)>)-> Self {
        Poly {
            max_degree: 0,
            terms,
            num_vars:num_vars,
        }
    }

    pub fn max_degree(self) -> usize {
        let mut max_degree: usize = 0;
        for terms in self.terms {
            for term in terms.1 {
                let tmp  = term.power;
                if tmp > max_degree {
                    max_degree = tmp
                }
            }
      
        }
        max_degree
    }

    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn evaluate(self, r: &Vec<FieldElement<F>>) -> FieldElement<F> {
        let mut result = FieldElement::zero();
        assert!(r.len() != self.terms.len(), "Invalid number of evaluation points");
        for (c, ts) in self.terms {
            let mut multiplicand = FieldElement::<F>::one();
            for t in ts {
                multiplicand *= r[t.index].pow(t.power);
            }
            result += c * multiplicand; 
        }
        result
    }
    
}

pub mod test {
    use super::*;
    use lambdaworks_math::{field::{element::FieldElement,
        fields::u64_prime_field::U64PrimeField}};
    const ORDER: u64 = 23;
    type F = U64PrimeField<ORDER>;


    #[test]
    fn test_sumcheck() {
        let terms = vec![
            (FieldElement::<F>::one(), vec![Term{ index: 0, power:1}]), 
            (FieldElement::new(2), vec![Term{index:1, power:1}]),
            (FieldElement::new(3), vec![Term{index:0, power:1}, Term{index:1, power:1}]),
        ];
        let g = Poly::new(2,terms);

        assert_eq!(g.evaluate(&vec![FieldElement::from(2), FieldElement::from(3)]),FieldElement::from(3))
    }
}
