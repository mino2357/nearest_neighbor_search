use apng::{load_dynamic_image, Encoder, Frame, PNGImage};
use plotters::prelude::*;
use rand::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::Path;

mod gen_grid2d;
mod kd_tree;

#[allow(dead_code)]
fn draw_graph(i: usize, vec: &gen_grid2d::Points2D, boundary: &gen_grid2d::Points2D) {
    let out_file_name = format!("{:04}", i).to_string() + ".png";

    let root = BitMapBackend::new(&out_file_name, (1080, 1080)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        //.caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-1.05..1.05, -1.05..1.05)
        .unwrap();

    // chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(PointSeries::of_element(
            (0..vec.points.len()).map(|i| (vec.points[i].x, vec.points[i].y)),
            2,
            ShapeStyle::from(&RED).filled(),
            &|coord, size, style| EmptyElement::at(coord) + Circle::new((0, 0), size, style),
        ))
        .unwrap();

    chart
        .draw_series(PointSeries::of_element(
            (0..boundary.points.len()).map(|i| (boundary.points[i].x, boundary.points[i].y)),
            2,
            ShapeStyle::from(&BLUE).filled(),
            &|coord, size, style| EmptyElement::at(coord) + Circle::new((0, 0), size, style),
        ))
        .unwrap();

    root.present().unwrap();
}

#[allow(dead_code)]
fn gen_apng(num: usize) {
    let mut files = vec![];

    for i in 0..num {
        files.push(format!("{:04}", i).to_string() + ".png");
    }

    let mut png_images: Vec<PNGImage> = Vec::new();

    for f in files.iter() {
        let mut file = File::open(f).unwrap();
        let mut buffer = vec![];
        file.read_to_end(&mut buffer).unwrap();
        let img = image::load_from_memory(&buffer).unwrap();
        png_images.push(load_dynamic_image(img).unwrap());
    }

    let path = Path::new(r"graph.png");
    let mut out = BufWriter::new(File::create(path).unwrap());

    let config = apng::create_config(&png_images, None).unwrap();
    let mut encoder = Encoder::new(&mut out, config).unwrap();

    for image in png_images.iter() {
        let frame = Frame {
            delay_num: Some(1),
            delay_den: Some(10),
            ..Default::default()
        };
        encoder.write_frame(image, frame).unwrap();
    }

    match encoder.finish_encode() {
        Ok(_n) => println!("success"),
        Err(err) => eprintln!("{}", err),
    }
}

fn main() {
    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let num_point: usize = 400;
    let num_boundary: usize = (std::f64::consts::PI * (num_point as f64).powf(0.5)) as usize;

    let mut vec = gen_grid2d::Points2D::new();
    let mut boundary = gen_grid2d::Points2D::new();

    loop {
        let x_r = 2.0 * (rng.gen::<f64>() - 0.5);
        let y_r = 2.0 * (rng.gen::<f64>() - 0.5);
        if x_r * x_r + y_r * y_r < 1.0 {
            vec.push(x_r, y_r);
        }
        if vec.points.len() == num_point {
            break;
        }
    }

    let radius = 4.0 * std::f64::consts::PI / (num_boundary as f64);
    println!("radius = {}", radius);

    let max_counter = 50;

    let mut tree = kd_tree::KDTree::construct_kd_tree(&vec);

    for i in 0..num_boundary {
        let t = 2.0 * std::f64::consts::PI * i as f64 / num_boundary as f64;
        boundary.push(t.cos(), t.sin());
    }

    for i in 0..max_counter {
        draw_graph(i, &vec, &boundary);
        for _ in 0..100 {
            vec.euler_step_by_near_points(&boundary, &tree, radius);
            //vec.euler_step(&boundary);
            tree = kd_tree::KDTree::construct_kd_tree(&vec);
        }
        println!("{} / {}", i, max_counter - 1);
    }

    gen_apng(max_counter);
}
