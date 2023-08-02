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
    pub fn new(vector: &gen_grid2d::Grid2D, id_: usize) -> Self {
        Self {
            id: id_,
            position: vector.clone(),
            left: None,
            right: None,
        }
    }

    #[allow(dead_code)]
    fn insert(&mut self, point: &gen_grid2d::Grid2D, mut depth: i32, id: usize) {
        //println!("{:#?}", &point);

        let axis = depth % 2;

        match axis {
            0 => {
                if self.position.x > point.x {
                    if let Some(left_node) = &mut self.left {
                        depth += 1;
                        left_node.insert(point, depth, id);
                    } else {
                        let node = Self::new(point, id);
                        self.left = Some(Box::new(node));
                    }
                    /*
                    let left_node = Some(&mut self.left);
                    match left_node {
                        Some(_) => {
                            depth += 1;
                            self.insert(point, depth, id);
                        }
                        _ => {
                            let node = Self::new(point, id);
                            self.left = Some(Box::new(node));
                        }
                    }
                    */
                } else {
                    if let Some(right_node) = &mut self.right {
                        depth += 1;
                        right_node.insert(point, depth, id);
                    } else {
                        let node = Self::new(point, id);
                        self.right = Some(Box::new(node));
                    }
                    /*
                    let right_node = Some(&mut self.right);
                    match right_node {
                        Some(_) => {
                            depth += 1;
                            self.insert(point, depth, id);
                        }
                        _ => {
                            let node = Self::new(point, id);
                            self.right = Some(Box::new(node));
                        }
                    }
                    */
                }
            }
            1 => {
                if self.position.y > point.y {
                    /*let left_node = Some(&mut self.left);
                    match left_node {
                        Some(_) => {
                            depth += 1;
                            self.insert(point, depth, id);
                        }
                        None => {
                            let node = Self::new(point, id);
                            self.left = Some(Box::new(node));
                        }
                    }*/
                    if let Some(left_node) = &mut self.left {
                        depth += 1;
                        left_node.insert(point, depth, id);
                    } else {
                        let node = Self::new(point, id);
                        self.left = Some(Box::new(node));
                    }
                } else {
                    /*
                    let right_node = Some(&mut self.right);
                    match right_node {
                        Some(_) => {
                            depth += 1;
                            self.insert(point, depth, id);
                        }
                        None => {
                            let node = Self::new(point, id);
                            self.right = Some(Box::new(node));
                        }
                    }
                    */
                    if let Some(right_node) = &mut self.right {
                        depth += 1;
                        right_node.insert(point, depth, id);
                    } else {
                        let node = Self::new(point, id);
                        self.right = Some(Box::new(node));
                    }
                }
            }
            _ => {}
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
    pub fn construct_kd_tree(vec: &mut gen_grid2d::Points2D) -> KDTree {
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
        let num_point: usize = 600;

        let mut vec = gen_grid2d::Points2D::new();

        for _ in 0..num_point {
            let x_r = 2.0 * (rng.gen::<f64>() - 0.5);
            let y_r = 2.0 * (rng.gen::<f64>() - 0.5);
            vec.push(x_r, y_r);
        }

        let tree = KDTree::construct_kd_tree(&mut vec);

        //println!("{:#?}", &tree);
        println!("{}", &tree.depth());
        println!("{}", &tree.size());
        assert_eq!(0, 1);
    }
}