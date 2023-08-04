use super::gen_grid2d;

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct KDTree
{
    pub id: usize,
    pub position: gen_grid2d::Grid2D,
    pub left: Option<Box<KDTree>>,
    pub right: Option<Box<KDTree>>,
}

impl KDTree {
    #[allow(dead_code)]
    fn new(vector: &gen_grid2d::Grid2D, id_: usize) -> Self {
        Self {
            id: id_,
            position: vector.clone(),
            left: None,
            right: None,
        }
    }

    #[allow(dead_code)]
    pub fn is_leaf(&self) -> bool {
        let left = &self.left;
        let right = &self.right;
        match (left, right) {
            (None, None) => {
                true
            }
            _ => {
                false
            }
        }
    }

    #[allow(dead_code)]
    pub fn search_points_id(&self, x: &gen_grid2d::Grid2D, radius: f64, near: &mut Vec<usize>, mut depth: i32) {

        let axis = depth % 2;
        //println!("id: {}, depth: {}, position: {:?}, {}", self.id, depth, self.position ,self.position.distance_square(x).sqrt());
        //println!("{:?}, {:?}", x, self.position);
        let r_self = self.position.distance_square(x).sqrt();
        if r_self < radius {
            near.push(self.id);
        }

        match axis {
            0 => {
                if self.position.x < x.x - radius {
                    match &self.right {
                        Some(right_node) => {
                            depth += 1;
                            right_node.search_points_id(x, radius, near, depth);
                        }
                        None => {}
                    }
                } else if x.x + radius < self.position.x {
                    match &self.left {
                        Some(left_node) => {
                            depth += 1;
                            left_node.search_points_id(x, radius, near, depth);
                        }
                        None => {}
                    }
                } else {
                    match (&self.right, &self.left) {
                        (Some(left_node), Some(right_node)) => {
                            depth += 1;
                            left_node.search_points_id(x, radius, near, depth);
                            right_node.search_points_id(x, radius, near, depth);
                        }
                        (Some(left_node), None) => {
                            depth += 1;
                            left_node.search_points_id(x, radius, near, depth);
                        }
                        (None, Some(right_node)) => {
                            depth += 1;
                            right_node.search_points_id(x, radius, near, depth);
                        }
                        (None, None) => {}
                    }
                }
            }
            _ => {
                if self.position.y < x.y - radius {
                    match &self.right {
                        Some(right_node) => {
                            depth += 1;
                            right_node.search_points_id(x, radius, near, depth);
                        }
                        None => {}
                    }
                } else if x.y + radius < self.position.y {
                    match &self.left {
                        Some(left_node) => {
                            depth += 1;
                            left_node.search_points_id(x, radius, near, depth);
                        }
                        None => {}
                    }
                } else {
                    match (&self.right, &self.left) {
                        (Some(left_node), Some(right_node)) => {
                            depth += 1;
                            left_node.search_points_id(x, radius, near, depth);
                            right_node.search_points_id(x, radius, near, depth);
                        }
                        (Some(left_node), None) => {
                            depth += 1;
                            left_node.search_points_id(x, radius, near, depth);
                        }
                        (None, Some(right_node)) => {
                            depth += 1;
                            right_node.search_points_id(x, radius, near, depth);
                        }
                        (None, None) => {}
                    }
                }
            }
        }
    }

    #[allow(dead_code)]
    fn insert(&mut self, point: &gen_grid2d::Grid2D, mut depth: i32, id: usize) {
        //println!("depth: {:?}", depth);
        //println!("point: {:#?}", point);

        let axis = depth % 2;

        match axis {
            0 => {
                if self.position.x > point.x {
                    match &mut self.left {
                        Some(left_node) => {
                            depth += 1;
                            left_node.insert(point, depth, id);
                        }
                        None => {
                            let node = Self::new(point, id);
                            self.left = Some(Box::new(node));
                        }
                    }
                } else {
                    match &mut self.right {
                        Some(right_node) => {
                            depth += 1;
                            right_node.insert(point, depth, id);
                        }
                        None => {
                            let node = Self::new(point, id);
                            self.right = Some(Box::new(node));
                        }
                    }
                }
            }
            _ => {
                if self.position.y > point.y {
                    match &mut self.left {
                        Some(left_node) => {
                            depth += 1;
                            left_node.insert(point, depth, id);
                        }
                        None => {
                            let node = Self::new(point, id);
                            self.left = Some(Box::new(node));
                        }
                    }
                }else {
                    match &mut self.right {
                        Some(right_node) => {
                            depth += 1;
                            right_node.insert(point, depth, id);
                        }
                        None => {
                            let node = Self::new(point, id);
                            self.right = Some(Box::new(node));
                        }
                    }
                }
            }
        }
    }

    #[allow(dead_code)]
    fn create_kd_tree(&mut self, vec: &gen_grid2d::Points2D) -> Self {
        let depth = 0;
        for i in 1..vec.points.len() {
            self.insert(&vec.points[i], depth, i);
        }
        self.clone()
    }

    #[allow(dead_code)]
    pub fn depth(&self) -> i32 {
        match (&self.left, &self.right) {
            (Some(left), Some(right)) => {
                1 + (left.depth()).max(right.depth())
            }
            (None, None) => {
                0
            }
            (None, Some(right)) => {
                1 + right.depth()
            }
            (Some(left), None) => {
                1 + left.depth()
            }
        }
    }

    #[allow(dead_code)]
    pub fn construct_kd_tree(vec: &gen_grid2d::Points2D) -> KDTree {
        let mut tree = KDTree::new(&gen_grid2d::Grid2D::new(vec.points[0].x, vec.points[0].y), 0);
        tree.create_kd_tree(&vec)
    }

    #[allow(dead_code)]
    pub fn size(&self) -> i32 {
        match (&self.left, &self.right) {
            (Some(left), Some(right)) => {
                1 + left.size() + right.size()
            }
            (None, None) => {
                1
            }
            (None, Some(right)) => {
                1 + right.size()
            }
            (Some(left), None) => {
                1 + left.size()
            }
        }
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn tree_1() {
        use rand::prelude::*;
        let seed: [u8; 32] = [1; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        let num_point: usize = 600;

        let mut vec = gen_grid2d::Points2D::new();

        for _ in 0..num_point {
            let x_r = 2.0 * (rng.gen::<f64>() - 0.5);
            let y_r = 2.0 * (rng.gen::<f64>() - 0.5);
            vec.push(x_r, y_r);
        }

        let tree = KDTree::construct_kd_tree(&mut vec);

        assert_eq!(tree.depth(), 20);
        assert_eq!(tree.size(), 600);
    }

    #[test]
    fn tree_2() {
        use rand::prelude::*;
        let seed: [u8; 32] = [1; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        let num_point: usize = 10;

        let mut vec = gen_grid2d::Points2D::new();

        for _ in 0..num_point {
            let x_r = 2.0 * (rng.gen::<f64>() - 0.5);
            let y_r = 2.0 * (rng.gen::<f64>() - 0.5);
            vec.push(x_r, y_r);
        }

        let tree = KDTree::construct_kd_tree(&mut vec);

        //println!("kd-tree depth: {}", &tree.depth());
        //println!("kd-tree size: {}", &tree.size());
        let center = gen_grid2d::Grid2D{ x: 0.0, y: 0.0};
        let radius = 0.4;
        //println!("radius: {}", radius);
        let mut near = vec![0 as usize; 0];
        tree.search_points_id(&center, radius, &mut near, 0);

        assert_eq!(near, [1 as usize, 9 as usize].to_vec());

    }

    #[test]
    fn tree_3() {
        use rand::prelude::*;
        let seed: [u8; 32] = [1; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        let num_point: usize = 10;

        let mut vec = gen_grid2d::Points2D::new();

        for _ in 0..num_point {
            let x_r = 2.0 * (rng.gen::<f64>() - 0.5);
            let y_r = 2.0 * (rng.gen::<f64>() - 0.5);
            vec.push(x_r, y_r);
        }

        let tree = KDTree::construct_kd_tree(&mut vec);

        //println!("kd-tree depth: {}", &tree.depth());
        //println!("kd-tree size: {}", &tree.size());
        let center = gen_grid2d::Grid2D{ x: 0.4, y: 0.3};
        let radius = 0.5;
        //println!("radius: {}", radius);
        let mut near = vec![0 as usize; 0];
        tree.search_points_id(&center, radius, &mut near, 0);

        //println!("{:?}", near);

        assert_eq!(near, [1 as usize, 2 as usize, 6 as usize, 9 as usize, 5 as usize].to_vec());

    }
}