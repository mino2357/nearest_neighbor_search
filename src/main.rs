use apng::{load_dynamic_image, Encoder, Frame, PNGImage};
use plotters::prelude::*;
use rand::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::Path;

#[allow(dead_code)]
#[derive(Debug, Clone)]
struct Grid2D {
    x: f64,
    y: f64,
}

#[derive(Debug)]
struct Points2D {
    pub points: Vec<Grid2D>,
}

impl Points2D {
    pub fn new() -> Self {
        Points2D { points: vec![] }
    }

    fn push(&mut self, x_r: f64, y_r: f64) {
        self.points.push(Grid2D { x: x_r, y: y_r });
    }

    fn euler_step(&mut self, boundary: &Points2D) {
        let dt = 1.0e-3;
        for i in 0..self.points.len() {
            //println!("{}", self.points[i].x);
            let x = self.points[i].x - dt * self.lennard_jones_potential_deriv(i, boundary).x;
            let y = self.points[i].y - dt * self.lennard_jones_potential_deriv(i, boundary).y;
            //let x = self.points[i].x - dt * self.simple_potential_deriv(i, &boundary).x;
            //let y = self.points[i].y - dt * self.simple_potential_deriv(i, &boundary).y;
            if x * x + y * y < 1.0 {
                self.points[i].x = x;
                self.points[i].y = y;
            } else {
                self.points[i].x = -0.5 * x;
                self.points[i].y = -0.5 * y;
            }
            //println!("{}", self.points[i].x);
        }
    }

    #[allow(dead_code)]
    fn lennard_jones_potential_deriv(&mut self, index: usize, boundary: &Points2D) -> Grid2D {
        let epsilon: f64 = 1.0;
        let sigma: f64 = 0.9 * (2.0_f64).powf(-1.0 / 6.0) / (self.points.len() as f64).powf(0.5);
        let mut f_x = 0.0;
        let mut f_y = 0.0;
        for i in 0..self.points.len() {
            if i != index {
                let dx = self.points[index].x - self.points[i].x;
                let dy = self.points[index].y - self.points[i].y;
                let r = (dx.powf(2.0) + dy.powf(2.0)).powf(0.5);
                f_x += 4.0
                    * epsilon
                    * (12.0 * sigma.powf(12.0) / r.powf(13.0)
                        - 6.0 * sigma.powf(6.0) / r.powf(7.0))
                    * dx;
                f_y += 4.0
                    * epsilon
                    * (12.0 * sigma.powf(12.0) / r.powf(13.0)
                        - 6.0 * sigma.powf(6.0) / r.powf(7.0))
                    * dy;
            }
        }
        for i in 0..boundary.points.len() {
            let dx = self.points[index].x - boundary.points[i].x;
            let dy = self.points[index].y - boundary.points[i].y;
            let r = (dx.powf(2.0) + dy.powf(2.0)).powf(0.5);
            f_x += 4.0
                * epsilon
                * (12.0 * sigma.powf(12.0) / r.powf(13.0) - 6.0 * sigma.powf(6.0) / r.powf(7.0))
                * dx;
            f_y += 4.0
                * epsilon
                * (12.0 * sigma.powf(12.0) / r.powf(13.0) - 6.0 * sigma.powf(6.0) / r.powf(7.0))
                * dy;
        }
        Grid2D { x: f_x, y: f_y }
    }

    #[allow(dead_code)]
    fn simple_potential_deriv(&mut self, index: usize, boundary: &Points2D) -> Grid2D {
        let a = 1.0e-5;
        let mut f_x = 0.0;
        let mut f_y = 0.0;
        for i in 0..self.points.len() {
            if i != index {
                let dx = self.points[index].x - self.points[i].x;
                let dy = self.points[index].y - self.points[i].y;
                let r = (dx.powf(2.0) + dy.powf(2.0)).powf(0.5);
                f_x -= a * 1.0 / r.powf(2.0) * dx;
                f_y -= a * 1.0 / r.powf(2.0) * dy;
            }
        }
        for i in 0..boundary.points.len() {
            let dx = self.points[index].x - boundary.points[i].x;
            let dy = self.points[index].y - boundary.points[i].y;
            let r = (dx.powf(2.0) + dy.powf(2.0)).powf(0.5);
            f_x -= a * 1.0 / r.powf(2.0) * dx;
            f_y -= a * 1.0 / r.powf(2.0) * dx;
        }
        Grid2D { x: f_x, y: f_y }
    }
}

fn draw_graph(i: usize, vec: &Points2D, boundary: &Points2D) {
    let out_file_name = format!("{:04}", i).to_string() + ".png";

    let root = BitMapBackend::new(&out_file_name, (1440, 1440)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        //.caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-1.05..1.05, -1.05..1.05)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

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

fn main() {
    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let num_point: usize = 1000;
    let num_boundary = (std::f64::consts::PI * (num_point as f64).powf(0.5)) as usize;

    let mut vec = Points2D::new();
    let mut boundary = Points2D::new();

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

    for i in 0..num_boundary {
        let t = 2.0 * std::f64::consts::PI * i as f64 / num_boundary as f64;
        boundary.push(t.cos(), t.sin());
    }

    let max_counter = 125;

    for i in 0..max_counter {
        //println!("{:?}", vec.points);
        //println!("{:?}", vec.lennard_jones_potential_deriv(0, &boundary));
        draw_graph(i, &vec, &boundary);
        for _ in 0..10 {
            vec.euler_step(&boundary);
        }
        println!("{}", i);
    }
    let mut files = vec![];

    for i in 0..max_counter {
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
            delay_den: Some(30), // 2, 3, 4, 5, 6, 7
            ..Default::default()
        };
        encoder.write_frame(image, frame).unwrap();
    }

    match encoder.finish_encode() {
        Ok(_n) => println!("success"),
        Err(err) => eprintln!("{}", err),
    }
}
