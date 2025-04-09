use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyfunction]
fn n_body_ode_system(
    time: f64,
    y: Vec<f64>,
    masses: Vec<f64>,
    g: f64,
) -> Vec<f64> {
    let num = masses.len();
    let mut bodies = Vec::with_capacity(num);

    // Extract position and velocity vectors
    for nth in 0..num {
        let vel = y[nth * 3..nth * 3 + 3].to_vec();
        let pos = y[num * 3 + nth * 3..num * 3 + nth * 3 + 3].to_vec();
        bodies.push((pos, vel));
    }

    let mut distances = vec![vec![0.0; num]; num];
    let mut accelerations = Vec::with_capacity(3 * num);
    let mut new_positions: Vec<f64> = Vec::with_capacity(3 * num);

    for first in 0..num {
        let mut acc = vec![0.0; 3];

        for second in 0..num {
            if first != second {
                if distances[first][second] == 0.0 {
                    let dist_squared: f64 = (0..3)
                        .map(|i| {
                            let d = bodies[second].0[i] - bodies[first].0[i];
                            d * d
                        })
                        .sum();

                    let dist = dist_squared.sqrt();
                    distances[first][second] = dist;
                    distances[second][first] = dist; // symmetric
                }

                let dist = distances[first][second];
                for i in 0..3 {
                    acc[i] += g * masses[second]
                        * (bodies[second].0[i] - bodies[first].0[i])
                        / dist.powi(3);
                }
            }
        }

        accelerations.extend(acc);
    }

    for nth in 0..num {
        new_positions.extend(&bodies[nth].1);
    }

    accelerations.extend(new_positions);
    accelerations
}

use pyo3::prelude::*;

#[pyfunction]
pub fn rkf45(
    y0: Vec<f64>,
    masses: Vec<f64>,
    t0: f64,
    t_end: f64,
    mut h: f64,
    g: f64,
    acceptable_error: f64,
    upper_bound_fraction: f64,
    lower_bound_fraction: f64,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let mut t_values = vec![];
    let mut y_values = vec![];

    let original_step_size = h;
    let mut t = t0;
    let mut y = y0.clone();
    let mut record_previous_data = true;
    let mut old_h;

    let a = [0.0, 0.0, 2.0 / 9.0, 1.0 / 3.0, 3.0 / 4.0, 1.0, 5.0 / 6.0];

    let b = [
        [0.0; 6],
        [0.0; 6],
        [0.0, 2.0 / 9.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0 / 12.0, 1.0 / 4.0, 0.0, 0.0, 0.0],
        [0.0, 69.0 / 128.0, -243.0 / 128.0, 135.0 / 64.0, 0.0, 0.0],
        [0.0, -17.0 / 12.0, 27.0 / 4.0, -27.0 / 5.0, 16.0 / 15.0, 0.0],
        [0.0, 65.0 / 432.0, -5.0 / 16.0, 13.0 / 16.0, 4.0 / 27.0, 5.0 / 144.0],
    ];

    let ch = [0.0, 47.0 / 450.0, 0.0, 12.0 / 25.0, 32.0 / 225.0, 1.0 / 30.0, 6.0 / 25.0];
    let ct = [1.0 / 150.0, 0.0, -3.0 / 100.0, 16.0 / 75.0, 1.0 / 20.0, -6.0 / 25.0];

    while t <= t_end {
        if record_previous_data {
            t_values.push(t);
            y_values.push(y.clone());
        }
        record_previous_data = true;

        let mut k = vec![vec![0.0; y.len()]; 6];
        let mut temp = vec![0.0; y.len()];

        // k1
        let mut f_t_y = n_body_ode_system(t, y.clone(), masses.clone(), g);
        for i in 0..y.len() {
            k[0][i] = h * f_t_y[i];
        }

        // k2 to k6
        for stage in 1..6 {
            for i in 0..y.len() {
                temp[i] = y[i];
                for j in 0..stage {
                    temp[i] += b[stage + 1][j + 1] * k[j][i];
                }
            }
            f_t_y = n_body_ode_system(t + a[stage + 1] * h, temp.clone(), masses.clone(), g);
            for i in 0..y.len() {
                k[stage][i] = h * f_t_y[i];
            }
        }

        // Update y
        for i in 0..y.len() {
            for j in 0..6 {
                y[i] += ch[j + 1] * k[j][i];
            }
        }

        // Estimate truncation error
        let mut truncation_error = 0.0;
        for i in 0..y.len() {
            let mut trunc_temp = 0.0;
            for j in 0..6 {
                trunc_temp += ct[j] * k[j][i];
            }
            truncation_error += trunc_temp * trunc_temp;
        }
        truncation_error = truncation_error.sqrt();

        // Adjust step size
        old_h = h;
        if truncation_error > acceptable_error {
            record_previous_data = false;
        }
        h *= 0.9 * (acceptable_error / truncation_error).powf(0.2);

        if h < lower_bound_fraction * original_step_size {
            h = lower_bound_fraction * original_step_size;
            record_previous_data = true;
        } else if h > upper_bound_fraction * original_step_size {
            h = upper_bound_fraction * original_step_size;
            record_previous_data = true;
        }

        if record_previous_data {
            t += old_h;
        }
    }

    Ok((t_values, y_values))
}


#[pymodule]
fn rusty_orbital_dynamics(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(n_body_ode_system, m)?)?;
    m.add_function(wrap_pyfunction!(rkf45, m)?)?;
    Ok(())
}
