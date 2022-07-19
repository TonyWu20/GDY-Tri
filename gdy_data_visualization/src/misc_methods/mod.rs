pub fn meshgrid<T>(x: Vec<T>, y: Vec<T>) -> (Vec<T>, Vec<T>)
where
    T: Clone,
{
    let nx = x.len();
    let ny = y.len();
    let yv: Vec<T> = (0..nx).into_iter().flat_map(|_| y.to_vec()).collect();
    let xv: Vec<T> = x.into_iter().flat_map(|value| vec![value; ny]).collect();
    (xv, yv)
}
