use std::collections::HashSet;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::Deserialize;
use serde_pyobject::from_pyobject;

mod buildup;
mod csg;
mod core;
mod source;

#[derive(Debug, Deserialize)]
struct CellInput {
    name: String,
    material_index: u32,
    csg_operations: Vec<i64>,
}

#[derive(Debug, Deserialize)]
struct PrimitveInput {
    shape: String,
    data: Vec<f64>
}

#[derive(Debug, Deserialize)]
struct MaterialInput {
    name: String,
    linear_attenuation_coefficient: Vec<f64>
}

#[derive(Debug, Deserialize, Clone)]
struct SourceInput {
    factor_by_group: Vec<f64>,
    source_type: String,
    fdata: Vec<f64>,
    udata: Vec<usize>,
}

#[derive(Debug, Deserialize)]
struct ToyPkInput {
    primitives: Vec<PrimitveInput>,
    cells: Vec<CellInput>,
    materials: Vec<MaterialInput>,
    source: SourceInput,
    tally_points: Vec<Vec<f64>>,
}


#[pyfunction]
fn run(input_dict: &PyDict) -> PyResult<()> {
    let toypk_input: ToyPkInput = from_pyobject(input_dict).unwrap();
    // println!("{:?}", toypk_input);
    let g = core::Geometry {
        primitives: toypk_input.primitives.iter().filter_map(|pi| 
            match &pi.shape[..] {
                "Sphere" => Some(Box::new(csg::Sph::from(&pi.data)) as Box<dyn csg::Primitive>),
                "Rpp" => Some(Box::new(csg::Rpp::from(&pi.data)) as Box<dyn csg::Primitive>),
                "Rcc" => Some(Box::new(csg::Rcc::from(&pi.data)) as Box<dyn csg::Primitive>),
                _ => None,
            }
        ).collect(),
        cells: toypk_input.cells.iter().map(|c|
            core::Cell {name: String::from(&c.name), material_index: c.material_index as usize, csg_operations: HashSet::from_iter(c.csg_operations.iter().cloned())}
        ).collect(),
        materials: toypk_input.materials.iter().map(|m|
            core::Material {name: String::from(&m.name), linear_attenuation_coefficient: m.linear_attenuation_coefficient.clone()}
        ).collect()
    };
    // TODO
    // implement buildup factor input conversion
    let bf = 2.0;
    let bu = buildup::BuildUpFactorUsed {
        d: vec![
            Box::new(buildup::TestingForm {bf}),
            Box::new(buildup::TestingForm {bf}),
        ]
    };
    // TODO
    // add new source type
    let s = toypk_input.source;
    let source = match &s.source_type[..] {
        "test" => source::TestSource::new(s.fdata[0], s.fdata[1], s.fdata[2], s.fdata[3], s.fdata[4], s.fdata[5], s.udata[0], s.udata[1], s.udata[2]),
        _ => panic!("invalid source type")
    };
    // TODO
    // calculate in parallel
    for p in toypk_input.tally_points {
        let p = csg::Point::new(p[0], p[1], p[2]);
        let mut cloned_source = source.clone();
        
        let result = core::run(&g, &mut cloned_source, &p, &s.factor_by_group, &bu);
        println!("{:?}", result);
    }
    Ok(())
}

#[pymodule]
fn toypk(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run, m)?)?;
    Ok(())
}