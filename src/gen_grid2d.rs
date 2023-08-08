use super::kd_tree;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Grid2D {
    pub x: f64,
    pub y: f64,
}

impl Grid2D {
    #[allow(dead_code)]
    pub fn new(x_: f64, y_: f64) -> Self {
        Grid2D { x: x_, y: y_ }
    }

    #[allow(dead_code)]
    pub fn distance_square(&self, point: &Grid2D) -> f64 {
        let dx2 = (self.x - point.x) * (self.x - point.x);
        let dy2 = (self.y - point.y) * (self.y - point.y);
        dx2 + dy2
    }
}

#[derive(Debug, Clone)]
pub struct Points2D {
    pub points: Vec<Grid2D>,
}

impl Points2D {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Points2D { points: vec![] }
    }

    #[allow(dead_code)]
    pub fn push(&mut self, x_r: f64, y_r: f64) {
        self.points.push(Grid2D { x: x_r, y: y_r });
    }

    #[allow(dead_code)]
    pub fn euler_step_by_near_points(
        &mut self,
        boundary: &Points2D,
        tree: &kd_tree::KDTree,
        radius: f64,
    ) {
        let dt = 1.0e-3;
        for i in 0..self.points.len() {
            let mut x = self.points[i].x
                - dt * self
                    .lennard_jones_potential_deriv_by_near_points(i, boundary, tree, radius)
                    .x;
            let mut y = self.points[i].y
                - dt * self
                    .lennard_jones_potential_deriv_by_near_points(i, boundary, tree, radius)
                    .y;
            if x * x + y * y < 1.0 {
                self.points[i].x = x;
                self.points[i].y = y;
            } else {
                loop {
                    let tmp = x;
                    x = -0.5 * y;
                    y = 0.5 * tmp;
                    if x * x + y * y < 1.0 {
                        self.points[i].x = x;
                        self.points[i].y = y;
                        break;
                    }
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn euler_step(&mut self, boundary: &Points2D) {
        let dt = 1.0e-3;
        for i in 0..self.points.len() {
            let mut x = self.points[i].x - dt * self.lennard_jones_potential_deriv(i, boundary).x;
            let mut y = self.points[i].y - dt * self.lennard_jones_potential_deriv(i, boundary).y;
            if x * x + y * y < 1.0 {
                self.points[i].x = x;
                self.points[i].y = y;
            } else {
                loop {
                    let tmp = x;
                    x = -0.5 * y;
                    y = 0.5 * tmp;
                    if x * x + y * y < 1.0 {
                        self.points[i].x = x;
                        self.points[i].y = y;
                        break;
                    }
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn lennard_jones_potential_deriv_by_near_points(
        &mut self,
        index: usize,
        boundary: &Points2D,
        tree: &kd_tree::KDTree,
        radius: f64,
    ) -> Grid2D {
        let epsilon: f64 = 1.0;
        let sigma: f64 = 0.9 * (2.0_f64).powf(-1.0 / 6.0) / (self.points.len() as f64).sqrt();
        let mut f_x = 0.0;
        let mut f_y = 0.0;
        let sigma2 = sigma * sigma;
        let sigma4 = sigma2 * sigma2;
        let sigma8 = sigma4 * sigma4;
        let sigma6 = sigma2 * sigma4;
        let sigma12 = sigma8 * sigma4;
        let near = tree.neighbor_search(&self.points[index], radius);
        for k in near.iter() {
            if index != *k {
                let dx = self.points[index].x - self.points[*k].x;
                let dy = self.points[index].y - self.points[*k].y;
                let r = (dx * dx + dy * dy).sqrt();
                let r2 = r * r;
                let r4 = r2 * r2;
                let r8 = r4 * r4;
                let r13 = r8 * r4 * r;
                let r7 = r2 * r4 * r;
                f_x += 4.0
                    * epsilon
                    * (12.0 * sigma12 / r13
                        - 6.0 * sigma6 / r7)
                    * dx;
                f_y += 4.0
                    * epsilon
                    * (12.0 * sigma12 / r13
                        - 6.0 * sigma6 / r7)
                    * dy;
            }
        }
        for i in 0..boundary.points.len() {
            let dx = self.points[index].x - boundary.points[i].x;
            let dy = self.points[index].y - boundary.points[i].y;
            let r = (dx * dx + dy * dy).sqrt();
            let r2 = r * r;
            let r4 = r2 * r2;
            let r8 = r4 * r4;
            let r13 = r8 * r4 * r;
            let r7 = r2 * r4 * r;
            f_x += 4.0
                * epsilon
                * (12.0 * sigma12 / r13 - 6.0 * sigma6 / r7)
                * dx;
            f_y += 4.0
                * epsilon
                * (12.0 * sigma12 / r13 - 6.0 * sigma6 / r7)
                * dy;
        }
        Grid2D { x: f_x, y: f_y }
    }

    #[allow(dead_code)]
    pub fn lennard_jones_potential_deriv(&mut self, index: usize, boundary: &Points2D) -> Grid2D {
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
}
