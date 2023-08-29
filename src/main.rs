use ndarray::Array2;
use ndarray_npy::{ReadNpyError, ReadNpyExt, WriteNpyExt};
use ndarray_stats::{QuantileExt,errors::MinMaxError};
use std::fs::File;
use std::collections::VecDeque;
use std::time::Instant;
use std::env;
use std::io::BufWriter;

#[derive(Debug,Copy,Clone)]
struct Cell {
    row: u32,
    col: u32,
    elevation: f64,
    is_boundary: bool,
}

fn planchon_darboux_alg(filename: &String,epsilon:f64,outputfile: &String) -> Result<(), ReadNpyError> {
    let reader = File::open(&filename)?;
    let dem:Array2<f64>= Array2::<f64>::read_npy(reader)?;
    let mut cells: Vec<Vec<Cell>> = Vec::new();
    let rows: usize = dem.shape()[0];
    let cols: usize = dem.shape()[1];
    
    //create the cell objects
    for row in 0..rows as usize{
        let mut cellrow: Vec<Cell> = Vec::new();
        for col in 0..cols as usize{
            let data = dem[[row,col]];
            let mut boundary_check: bool = false;
            if (row==0) || (col == 0) || (row == rows - 1) || (col == cols - 1) {
                boundary_check = true;
            }


            let cell = Cell{
                row:row as u32,
                col:col as u32,
                elevation:data as f64,
                is_boundary: boundary_check};
            cellrow.push(cell);
        }
        cells.push(cellrow);
    }

    let result: Result <&f64, MinMaxError>  = dem.max();
    let mut max_value : f64 = 3000.0;
    match result {
        Ok(value) => {max_value = *value},
        Err(err) => {println!("{}",err)},
    }

    let mut filled_dem: Array2<f64> = Array2::from_elem((rows,cols),max_value);
    
    let mut priority: VecDeque<Cell> = VecDeque::new();
    let mut queue: VecDeque<Cell> = VecDeque::new();
    fn dry_cell(cell:Cell, filled_dem: &mut Array2<f64>, priority: &mut VecDeque<Cell>){
        filled_dem[[cell.row as usize,cell.col as usize]] = cell.elevation;
        priority.push_back(cell);
    }
    
    //initialize filled_dem
    for row in 0..rows as usize {
        for col in 0..cols as usize {
            let cell: &Cell = &cells[row][col];
            if cell.is_boundary {
                dry_cell(*cell, &mut filled_dem, &mut priority);
            }
        }
    }

    println!("{},{}",filled_dem[[0,0]],cells[0][0].elevation);



    let offsets:[(i32,i32); 8] = [(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1),(0,-1)];
    let rowoffs: Vec<i32>= offsets.iter().map(|item| item.0).collect();
    let coloffs: Vec<i32>= offsets.iter().map(|item| item.1).collect();

    let mut iterations: u128 = 0;
    
    while (priority.len() > 0) || (queue.len() > 0) {
        if iterations % 1000 == 0 {
            println!("Iterations: {}, len_P: {}, len_Q: {}",iterations,priority.len(),queue.len());
        }

        let cell: Cell;
        if priority.len() > 0 {
            cell =  priority.pop_front().unwrap();
        } else {
            
            cell = queue.pop_front().unwrap();
            let cell_dem:f64 = dem[[cell.row as usize,cell.col as usize]];
            let cell_water:f64 = filled_dem[[cell.row as usize,cell.col as usize]];

            if cell_dem == cell_water {
                continue;
            }
        }

        let cell_water: f64 = filled_dem[[cell.row as usize, cell.col as usize]];

        for i in 0..8 as usize {
            let nextrow: i32 = (cell.row as i32) + rowoffs[i];
            let nextcol: i32 = (cell.col as i32) + coloffs[i];

            if (nextrow >= 0) && (nextrow < (rows as i32)- 1) && (nextcol >= 0) && (nextcol < (cols as i32) - 1){
                let neighbor: Cell = cells[nextrow as usize][nextcol as usize];
                let neighbor_elevation: f64 = neighbor.elevation;
                let neighbor_water: f64 = filled_dem[[neighbor.row as usize,neighbor.col as usize]];

                if neighbor_elevation >= neighbor_water{
                    continue;
                }

                if neighbor_elevation >= cell_water + epsilon {
                    dry_cell(neighbor,&mut filled_dem, &mut priority);
                } else if neighbor_water > cell_water + epsilon {
                    filled_dem[[neighbor.row as usize, neighbor.col as usize]] = cell_water + epsilon;
                    queue.push_back(neighbor);
                }


            }

        }

        
        iterations += 1;
    }

    println!("\n--------------------------------------------");
    println!("Sink-filling done.\nNumber of Iterations: {}\nepsilon: {}",iterations,epsilon);

    let writer = BufWriter::new(File::create(&outputfile)?);
    match filled_dem.write_npy(writer){
        Ok(_) => {println!("{} was written.",outputfile)},
        Err(err) => {println!("An error occured while writing {}. See below:\n{}",outputfile,err)},
    };

    Ok(())
}

fn main() {
    
    //command line arguments
    let args: Vec<String> = env::args().collect();
    
    //check if number of inputs is correct
    if args.len() != 4 {
        eprintln!("Usage: {} <filename.npy> <float> <outputfile.npy>", args[0]);
        std::process::exit(1);
    }

    //get filename
    let filename: &String = &args[1];

    //get output filename
    let outputfile: &String = &args[3];

    //get float argument and parse
    let str_epsilon: &String = &args[2];
    let parsed_epsilon: Result<f64, _> = str_epsilon.parse();

    //check if parsing was successful
    let epsilon = match parsed_epsilon {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Invalid epsilon value: {}",str_epsilon);
            std::process::exit(1);
        }
    };

    let start = Instant::now();
    planchon_darboux_alg(filename,epsilon,outputfile).expect("Failed to start sink-filling routine.\n");
    let elapsed = start.elapsed();

    println!("Elapsed time: {:.2?}",elapsed);
}