/**
Return meshgrid `XV` and `YV` from input Vec<T> of `x` and `y`
`XV` and `YV` are inherently of same dimension (nx, ny) and flattened before output.
# Arguments:
- x: `Vec<T>` - vector of data for x-axis values.
- y: `Vec<T>` - vector of data for y-axis values.
# Returns:
- (xv, yv) : `(Vec<T>, Vec<T>)` - xv is x repeated by ny times, yv is y[i] repeated with nx times
*/
pub fn meshgrid<T>(x: Vec<T>, y: Vec<T>) -> (Vec<T>, Vec<T>)
where
    T: Clone,
{
    let nx = x.len();
    let ny = y.len();
    let xv: Vec<T> = (0..ny).into_iter().flat_map(|_| x.to_vec()).collect();
    let yv: Vec<T> = y.into_iter().flat_map(|value| vec![value; nx]).collect();
    (xv, yv)
}
