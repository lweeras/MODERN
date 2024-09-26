# If desired, used setwd() to declare directory where data file should be stored
# setwd()

# Set the seed for data generation and declare the size and number of species to include in the data set
seed <- 2
num_row_col <- 5
num_species <- 3
species_probs <- c(0.5,0.33,0.1)

# Set the seed
set.seed(seed)

# Open a file to store the data
fileConn <- file("prob.dat", open = "w+")

# Initialize a counter for each species.
# Then randomly assign species to cells based on assigned probabilities and update the counters.
# Write assignments to file.
species_counts <- rep(0, num_species)
for (i in 1:num_row_col)
{
  for (j in 1:num_row_col)
  {
    for (k in 1:num_species)
    {
      r <- runif(1)
      if(r <= species_probs[k])
      {
        writeLines(paste(c(k-1, i-1, j-1), collapse = " ", sep = ""), fileConn)
        species_counts[k] <- species_counts[k] + 1
      }
    }
  }
}

# Write an empty line.
writeLines("", fileConn)

# Generate cost data.
# We assume that cost is impacted (somewhat) by costs of surrounding cells.
cost_mat <- matrix(rep(0, num_row_col^2), nrow = num_row_col)
order_mat <- matrix(rep(0, num_row_col^2), nrow = num_row_col)
unassigned_indices <- seq(0, num_row_col^2 - 1)
write_counter <- 0
while(length(unassigned_indices) > 0)
{
  # Find a cell that doesn't have an assigned cost
  index <- unassigned_indices[1]
  if(length(unassigned_indices) > 1)
  {
    index <- sample(unassigned_indices, 1)
  }
  unassigned_indices <- unassigned_indices[unassigned_indices != index]
  
  # Compute the average cost of adjacent cells that already have assigned costs.
  num_assigned_adjacent <- 0
  sum_adjacent <- 0
  if(ceiling((index+1)/num_row_col) > 1 && 
     cost_mat[ceiling((index+1)/num_row_col) - 1, index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1] > 0)
  {
    num_assigned_adjacent <- num_assigned_adjacent+1
    sum_adjacent <- sum_adjacent + cost_mat[ceiling((index+1)/num_row_col) - 1, index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1]
  }
  if(ceiling((index+1)/num_row_col) < num_row_col && 
     cost_mat[ceiling((index+1)/num_row_col) + 1, index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1] > 0)
  {
    num_assigned_adjacent <- num_assigned_adjacent+1
    sum_adjacent <- sum_adjacent + cost_mat[ceiling((index+1)/num_row_col) + 1, index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1]
  }
  if(index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1 > 1 && 
     cost_mat[ceiling((index+1)/num_row_col), index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1 - 1] > 0)
  {
    num_assigned_adjacent <- num_assigned_adjacent+1
    sum_adjacent <- sum_adjacent + cost_mat[ceiling((index+1)/num_row_col), index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1 - 1]
  }
  if(index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1 < num_row_col && 
     cost_mat[ceiling((index+1)/num_row_col), index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1 + 1] > 0)
  {
    num_assigned_adjacent <- num_assigned_adjacent+1
    sum_adjacent <- sum_adjacent + cost_mat[ceiling((index+1)/num_row_col), index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1 + 1]
  }
  mean_adjacent <- 0
  if(num_assigned_adjacent > 0)
  {
    mean_adjacent <- sum_adjacent/num_assigned_adjacent
  }
  
  # Assign a cost to the cell. 
  # If any adjacent cell has an assigned cost, the assigned cost should be "close" to the average of the assigned values.
  # Otherwise, use a triangular distribution from $0.1 to $1.9 with mode $1 (all in millions).
  if(mean_adjacent <= 0)
  {
    cost_mat[ceiling((index+1)/num_row_col), index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1] <- EnvStats::rtri(1,
                                                                                                                 min = 0.1,
                                                                                                                 max = 1.9,
                                                                                                                 mode = 1)
  }
  else
  {
    cost_mat[ceiling((index+1)/num_row_col), index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1] <- EnvStats::rtri(1,
                                                                                                                 min = max(c(.1, mean_adjacent - 0.5/num_assigned_adjacent)),
                                                                                                                 max = mean_adjacent + 0.5/num_assigned_adjacent,
                                                                                                                 mode = mean_adjacent)
  }
  order_mat[ceiling((index+1)/num_row_col), index - num_row_col*(ceiling((index+1)/num_row_col)-1) + 1] <- write_counter
  write_counter <- write_counter + 1
}

# Write cost assignments to file, but round to $0.25 increments.
for (i in 1:num_row_col)
{
  for (j in 1:num_row_col)
  {
    writeLines(paste(c(i-1, j-1, ceiling(cost_mat[i,j] / 0.25)*0.25), collapse = " ", sep = ""), fileConn)
  }
}

# Write a blank line.
writeLines("", fileConn)

# For each species, randomly generate the number cells that are required that contain that species.
# Write requirements to file.
species_req <- rep(0, num_species)
for(i in 1:num_species)
{
  species_req[i] <- round(runif(1, 1, min(c(species_counts[i],20))))
}
writeLines(paste(species_req, collapse = " ", sep = ""), fileConn)

# Write a blank line.
writeLines("", fileConn)

# Compute a rudimentary limit on total budget.
# Note that this will not always be feasible. In such cases, we manually increased the budget requirement to obtain feasibility. Specifically, we first round up to the next higher integer. If this is still not feasible, we increase by 1 until feasible.
avg_cost <- sum(cost_mat)/(num_row_col**2)
max_num_sites <- min(c(sum(species_req),num_row_col**2))
budget_limit <- round(runif(1, avg_cost*max(species_req), avg_cost*max_num_sites), 3)

writeLines(as.character(budget_limit), fileConn)

close(fileConn)