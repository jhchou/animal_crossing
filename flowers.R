library(tidyverse)

# Read flower genotype / phenotype / source data
# - data source, ACNH Flower Research: https://docs.google.com/spreadsheets/d/1rbYbQ0i3SuTu30KTma5dO4uuJWr_SjOZXA1l4UOIHWo
# - if data source file unavailable, use the code at the bottom to recreate data dataframe
data <- read.csv(file = 'flowers.csv', stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename_all(tolower) %>%
  mutate(
    source = ifelse(source == 1, 'seed', ifelse(source == 2, 'island', 'cross')),
    seedbag = (seedbag != 0),
    breedtrue = ((red != 1) & (yellow != 1) & (white != 1) & (brightness != 1)),
    genotype = paste0(red, yellow, white, brightness)
  ) %>%
  mutate_if(is.character, tolower) %>%
  mutate(genotype = ifelse(type == 'rose', genotype, substr(genotype, 1, nchar(genotype) - 1))) %>%
  mutate(genotype = paste(genotype, type)) %>%
  select(-id, -seedbag)

# Function to determine Mendelian results of genotype crosses
# - input: g1 / g2 genotypes: either 0, 1, or 2
# - return possible genotypes
cross_genotypes <- function(g1, g2) {
  g <- sort(c(g1, g2))
  if          (all(g == c(0,0))) { result <- c(0)
    } else if (all(g == c(0,1))) { result <- c(0, 1)
    } else if (all(g == c(0,2))) { result <- c(1)
    } else if (all(g == c(1,1))) { result <- c(0, 1, 1, 2)
    } else if (all(g == c(1,2))) { result <- c(1, 2)
    } else if (all(g == c(2,2))) { result <- c(2)
    }
  return(result)
}

# Results of cross including proportions
# - cross_flowers_single('tulip 001 020')
# - b1 and b2 = 0 for all except rose
cross_flowers_single <- function(input) {
  temp <- unlist(strsplit(input, ' '))
  type <- tolower(temp[1])
  genotype1 <- temp[2]
  genotype2 <- temp[3]
  
  g1 <- as.numeric(unlist(strsplit(genotype1, split = ''))) # vector of the red / yellow / white / brightness genotypes for 1st
  g2 <- as.numeric(unlist(strsplit(genotype2, split = ''))) # same for 2nd
  if((type == 'rose') & ((length(g1) != 4) | (length(g2) != 4))) { warning('Rose missing "brightness" genotype -- assuming 0', call. = FALSE) }
  
  r1 = g1[1]
  y1 = g1[2]
  w1 = g1[3]
  b1 = ifelse(length(g1) == 4, g1[4], 0) # if not a rose, brightness genotype is 0
  r2 = g2[1]
  y2 = g2[2]
  w2 = g2[3]
  b2 = ifelse(length(g2) == 4, g2[4], 0)
  
  result <- expand.grid(
    type       = type,
    red        = cross_genotypes(r1, r2),
    yellow     = cross_genotypes(y1, y2),
    white      = cross_genotypes(w1, w2),
    brightness = cross_genotypes(b1, b2),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    left_join(data %>% select(-source), by = c("type", "red", "yellow", "white", "brightness")) %>%
    mutate(total_n = n()) %>%
    group_by(genotype, type, red, yellow, white, brightness, color, breedtrue, total_n) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    group_by(color) %>%
    mutate(
      prop_cross = n/ total_n, # proportion this outcome, out of ALL possibilities of cross
      prop_color = n / sum(n)  # proportion this outcome, for this color phenotype
    ) %>%
    ungroup() %>%
    select(-n, -total_n) %>%
    arrange(color, desc(prop_color))
  
  return(result %>% select(genotype, color, breedtrue, prop_cross, prop_color))
}



cross_flowers <- function(input, single_cross = FALSE) {
  # Input any number of genotypes and generate results of all possible crosses
  # Combinations with repetitions: https://stackoverflow.com/questions/55738847/r-creating-combinations-with-replacement
  input <- unlist(strsplit(input, ' ')) # split string to vector of space-separated components
  type <- tolower(input[1]) # first element in flower type
  genotypes <- input[-1]    # remaining elements are genotypes
  
  if (length(genotypes) > 2) {
    if(single_cross) { warning('Too many genotypes -- displaying all possible crosses', call. = FALSE) }
    single_cross <- FALSE
  } else if (length(genotypes) == 1) {
    single_cross <- FALSE # grid expansion will fix this
  }
  
  if (single_cross) {
    result <- data.frame(list(Var1 = genotypes[1], Var2 = genotypes[2]))
  } else {
    result <- expand.grid(genotypes, genotypes, stringsAsFactors = FALSE) %>%
      mutate(key = paste(pmin(Var1, Var2), pmax(Var1, Var2), sep = "-")) %>%
      filter(!duplicated(key)) %>%
      select(-key)
  }
  
  result <- result %>%
    mutate(
      parent1 = paste(Var1, type),
      parent2 = paste(Var2, type)
    ) %>%
    left_join(data %>% select(genotype, color1 = color), by = c('parent1' = 'genotype')) %>% # color of first parent
    left_join(data %>% select(genotype, color2 = color), by = c('parent2' = 'genotype')) %>% # color of second parent
    mutate(
      cross = paste0(Var1, ' (', color1, ') x ', Var2, ' (', color2, ')')
    ) %>%
    mutate(input = paste(type, Var1, Var2)) %>%
    mutate(output = map(input, cross_flowers_single)) %>%
    unnest(output) %>%
    select(cross, genotype, color, breedtrue, prop_cross, prop_color)
  
  return(result)
}





# If 'flowers.csv' unavailable, the data can be reconstructed from the dput(data) as follows
# data <- structure(list(type = c("rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "rose", "rose", "rose", "rose", 
#   "rose", "rose", "rose", "rose", "tulip", "tulip", "tulip", "tulip", 
#   "tulip", "tulip", "tulip", "tulip", "tulip", "tulip", "tulip", 
#   "tulip", "tulip", "tulip", "tulip", "tulip", "tulip", "tulip", 
#   "tulip", "tulip", "tulip", "tulip", "tulip", "tulip", "tulip", 
#   "tulip", "tulip", "pansy", "pansy", "pansy", "pansy", "pansy", 
#   "pansy", "pansy", "pansy", "pansy", "pansy", "pansy", "pansy", 
#   "pansy", "pansy", "pansy", "pansy", "pansy", "pansy", "pansy", 
#   "pansy", "pansy", "pansy", "pansy", "pansy", "pansy", "pansy", 
#   "pansy", "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", 
#   "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", 
#   "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", 
#   "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", "cosmos", 
#   "lily", "lily", "lily", "lily", "lily", "lily", "lily", "lily", 
#   "lily", "lily", "lily", "lily", "lily", "lily", "lily", "lily", 
#   "lily", "lily", "lily", "lily", "lily", "lily", "lily", "lily", 
#   "lily", "lily", "lily", "hyacinth", "hyacinth", "hyacinth", "hyacinth", 
#   "hyacinth", "hyacinth", "hyacinth", "hyacinth", "hyacinth", "hyacinth", 
#   "hyacinth", "hyacinth", "hyacinth", "hyacinth", "hyacinth", "hyacinth", 
#   "hyacinth", "hyacinth", "hyacinth", "hyacinth", "hyacinth", "hyacinth", 
#   "hyacinth", "hyacinth", "hyacinth", "hyacinth", "hyacinth", "windflower", 
#   "windflower", "windflower", "windflower", "windflower", "windflower", 
#   "windflower", "windflower", "windflower", "windflower", "windflower", 
#   "windflower", "windflower", "windflower", "windflower", "windflower", 
#   "windflower", "windflower", "windflower", "windflower", "windflower", 
#   "windflower", "windflower", "windflower", "windflower", "windflower", 
#   "windflower", "mum", "mum", "mum", "mum", "mum", "mum", "mum", 
#   "mum", "mum", "mum", "mum", "mum", "mum", "mum", "mum", "mum", 
#   "mum", "mum", "mum", "mum", "mum", "mum", "mum", "mum", "mum", 
#   "mum", "mum"), red = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
#   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
#   2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
#   2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
#   2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 
#   1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
#   2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 
#   2L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 
#   1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
#   2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
#   2L, 2L, 2L, 2L), yellow = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
#   2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 
#   1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
#   2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 
#   2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 
#   2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 
#   1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 
#   0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 
#   0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 
#   2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 
#   1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 
#   1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 
#   0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 
#   2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 
#   2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 
#   1L, 1L, 2L, 2L, 2L), white = c(0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 
#   2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 
#   2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 
#   1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 
#   0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 
#   0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 
#   1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 
#   2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 
#   0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 
#   1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 
#   2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 
#   0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 
#   1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 
#   2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 
#   0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 
#   1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 
#   2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 
#   0L, 1L, 2L, 0L, 1L, 2L), brightness = c(0L, 1L, 2L, 0L, 1L, 2L, 
#   0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 
#   1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 
#   2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 
#   0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 
#   1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
#   0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), color = c("white", "white", 
#   "white", "white", "white", "white", "purple", "purple", "purple", 
#   "yellow", "yellow", "yellow", "white", "white", "white", "purple", 
#   "purple", "purple", "yellow", "yellow", "yellow", "yellow", "yellow", 
#   "yellow", "white", "white", "white", "red", "pink", "white", 
#   "red", "pink", "white", "red", "pink", "purple", "orange", "yellow", 
#   "yellow", "red", "pink", "white", "red", "pink", "purple", "orange", 
#   "yellow", "yellow", "orange", "yellow", "yellow", "red", "pink", 
#   "white", "black", "red", "pink", "black", "red", "pink", "black", 
#   "red", "pink", "orange", "orange", "yellow", "red", "red", "white", 
#   "black", "red", "purple", "orange", "orange", "yellow", "orange", 
#   "orange", "yellow", "blue", "red", "white", "white", "white", 
#   "white", "yellow", "yellow", "white", "yellow", "yellow", "yellow", 
#   "red", "pink", "white", "orange", "yellow", "yellow", "orange", 
#   "yellow", "yellow", "black", "red", "red", "black", "red", "red", 
#   "purple", "purple", "purple", "white", "white", "blue", "yellow", 
#   "yellow", "blue", "yellow", "yellow", "yellow", "red", "red", 
#   "blue", "orange", "orange", "orange", "yellow", "yellow", "yellow", 
#   "red", "red", "purple", "red", "red", "purple", "orange", "orange", 
#   "purple", "white", "white", "white", "yellow", "yellow", "white", 
#   "yellow", "yellow", "yellow", "pink", "pink", "pink", "orange", 
#   "orange", "pink", "orange", "orange", "orange", "red", "red", 
#   "red", "orange", "orange", "red", "black", "black", "red", "white", 
#   "white", "white", "yellow", "white", "white", "yellow", "yellow", 
#   "white", "red", "pink", "white", "orange", "yellow", "yellow", 
#   "orange", "yellow", "yellow", "black", "red", "pink", "black", 
#   "red", "pink", "orange", "orange", "white", "white", "white", 
#   "blue", "yellow", "yellow", "white", "yellow", "yellow", "yellow", 
#   "red", "pink", "white", "orange", "yellow", "yellow", "orange", 
#   "yellow", "yellow", "red", "red", "red", "blue", "red", "red", 
#   "purple", "purple", "purple", "white", "white", "blue", "orange", 
#   "orange", "blue", "orange", "orange", "orange", "red", "red", 
#   "blue", "pink", "pink", "pink", "orange", "orange", "orange", 
#   "red", "red", "purple", "red", "red", "purple", "pink", "pink", 
#   "purple", "white", "white", "purple", "yellow", "yellow", "white", 
#   "yellow", "yellow", "yellow", "pink", "pink", "pink", "yellow", 
#   "red", "pink", "purple", "purple", "purple", "red", "red", "red", 
#   "purple", "purple", "red", "green", "green", "red"), source = c("cross", 
#   "cross", "cross", "seed", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "seed", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "seed", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "island", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "island", "cross", 
#   "cross", "cross", "cross", "cross", "seed", "cross", "cross", 
#   "cross", "cross", "seed", "cross", "cross", "cross", "island", 
#   "cross", "cross", "cross", "cross", "island", "cross", "cross", 
#   "cross", "seed", "cross", "island", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "seed", "cross", "cross", "cross", 
#   "cross", "seed", "cross", "cross", "cross", "cross", "island", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "seed", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "island", 
#   "cross", "cross", "seed", "cross", "cross", "cross", "cross", 
#   "cross", "seed", "cross", "cross", "cross", "cross", "cross", 
#   "cross", "island", "cross", "cross", "cross", "seed", "cross", 
#   "cross", "cross", "island", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "seed", "cross", "cross", "cross", "seed", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "cross", "seed", "cross", 
#   "island", "cross", "island", "cross", "island", "cross", "cross", 
#   "seed", "cross", "cross", "cross", "cross", "seed", "cross", 
#   "cross", "cross", "island", "cross", "cross", "cross", "cross", 
#   "island", "cross", "cross", "cross", "seed", "cross", "island", 
#   "cross", "cross", "cross", "cross", "cross", "cross", "seed", 
#   "cross", "cross", "cross", "cross", "seed", "cross", "cross", 
#   "cross", "cross", "island", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "seed", "cross", "cross", "cross", "cross", 
#   "cross", "cross", "island", "cross", "cross", "seed", "cross", 
#   "cross", "cross", "cross", "seed", "cross", "cross", "cross", 
#   "cross", "cross", "cross", "cross", "island", "cross", "cross", 
#   "cross", "seed", "cross", "cross", "cross", "island", "cross", 
#   "cross", "cross", "cross"), breedtrue = c(TRUE, FALSE, TRUE, 
#   FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, 
#   FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, 
#   FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, 
#   FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, 
#   FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, 
#   FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, 
#   FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, 
#   FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, 
#   TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, 
#   TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, 
#   FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, 
#   TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
#   TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, 
#   FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, 
#   TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, 
#   FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, 
#   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, 
#   FALSE, FALSE, FALSE, TRUE, FALSE, TRUE), genotype = c("0000 rose", 
#   "0001 rose", "0002 rose", "0010 rose", "0011 rose", "0012 rose", 
#   "0020 rose", "0021 rose", "0022 rose", "0100 rose", "0101 rose", 
#   "0102 rose", "0110 rose", "0111 rose", "0112 rose", "0120 rose", 
#   "0121 rose", "0122 rose", "0200 rose", "0201 rose", "0202 rose", 
#   "0210 rose", "0211 rose", "0212 rose", "0220 rose", "0221 rose", 
#   "0222 rose", "1000 rose", "1001 rose", "1002 rose", "1010 rose", 
#   "1011 rose", "1012 rose", "1020 rose", "1021 rose", "1022 rose", 
#   "1100 rose", "1101 rose", "1102 rose", "1110 rose", "1111 rose", 
#   "1112 rose", "1120 rose", "1121 rose", "1122 rose", "1200 rose", 
#   "1201 rose", "1202 rose", "1210 rose", "1211 rose", "1212 rose", 
#   "1220 rose", "1221 rose", "1222 rose", "2000 rose", "2001 rose", 
#   "2002 rose", "2010 rose", "2011 rose", "2012 rose", "2020 rose", 
#   "2021 rose", "2022 rose", "2100 rose", "2101 rose", "2102 rose", 
#   "2110 rose", "2111 rose", "2112 rose", "2120 rose", "2121 rose", 
#   "2122 rose", "2200 rose", "2201 rose", "2202 rose", "2210 rose", 
#   "2211 rose", "2212 rose", "2220 rose", "2221 rose", "2222 rose", 
#   "000 tulip", "001 tulip", "002 tulip", "010 tulip", "011 tulip", 
#   "012 tulip", "020 tulip", "021 tulip", "022 tulip", "100 tulip", 
#   "101 tulip", "102 tulip", "110 tulip", "111 tulip", "112 tulip", 
#   "120 tulip", "121 tulip", "122 tulip", "200 tulip", "201 tulip", 
#   "202 tulip", "210 tulip", "211 tulip", "212 tulip", "220 tulip", 
#   "221 tulip", "222 tulip", "000 pansy", "001 pansy", "002 pansy", 
#   "010 pansy", "011 pansy", "012 pansy", "020 pansy", "021 pansy", 
#   "022 pansy", "100 pansy", "101 pansy", "102 pansy", "110 pansy", 
#   "111 pansy", "112 pansy", "120 pansy", "121 pansy", "122 pansy", 
#   "200 pansy", "201 pansy", "202 pansy", "210 pansy", "211 pansy", 
#   "212 pansy", "220 pansy", "221 pansy", "222 pansy", "000 cosmos", 
#   "001 cosmos", "002 cosmos", "010 cosmos", "011 cosmos", "012 cosmos", 
#   "020 cosmos", "021 cosmos", "022 cosmos", "100 cosmos", "101 cosmos", 
#   "102 cosmos", "110 cosmos", "111 cosmos", "112 cosmos", "120 cosmos", 
#   "121 cosmos", "122 cosmos", "200 cosmos", "201 cosmos", "202 cosmos", 
#   "210 cosmos", "211 cosmos", "212 cosmos", "220 cosmos", "221 cosmos", 
#   "222 cosmos", "000 lily", "001 lily", "002 lily", "010 lily", 
#   "011 lily", "012 lily", "020 lily", "021 lily", "022 lily", "100 lily", 
#   "101 lily", "102 lily", "110 lily", "111 lily", "112 lily", "120 lily", 
#   "121 lily", "122 lily", "200 lily", "201 lily", "202 lily", "210 lily", 
#   "211 lily", "212 lily", "220 lily", "221 lily", "222 lily", "000 hyacinth", 
#   "001 hyacinth", "002 hyacinth", "010 hyacinth", "011 hyacinth", 
#   "012 hyacinth", "020 hyacinth", "021 hyacinth", "022 hyacinth", 
#   "100 hyacinth", "101 hyacinth", "102 hyacinth", "110 hyacinth", 
#   "111 hyacinth", "112 hyacinth", "120 hyacinth", "121 hyacinth", 
#   "122 hyacinth", "200 hyacinth", "201 hyacinth", "202 hyacinth", 
#   "210 hyacinth", "211 hyacinth", "212 hyacinth", "220 hyacinth", 
#   "221 hyacinth", "222 hyacinth", "000 windflower", "001 windflower", 
#   "002 windflower", "010 windflower", "011 windflower", "012 windflower", 
#   "020 windflower", "021 windflower", "022 windflower", "100 windflower", 
#   "101 windflower", "102 windflower", "110 windflower", "111 windflower", 
#   "112 windflower", "120 windflower", "121 windflower", "122 windflower", 
#   "200 windflower", "201 windflower", "202 windflower", "210 windflower", 
#   "211 windflower", "212 windflower", "220 windflower", "221 windflower", 
#   "222 windflower", "000 mum", "001 mum", "002 mum", "010 mum", 
#   "011 mum", "012 mum", "020 mum", "021 mum", "022 mum", "100 mum", 
#   "101 mum", "102 mum", "110 mum", "111 mum", "112 mum", "120 mum", 
#   "121 mum", "122 mum", "200 mum", "201 mum", "202 mum", "210 mum", 
#   "211 mum", "212 mum", "220 mum", "221 mum", "222 mum")), class = c("tbl_df", 
#   "tbl", "data.frame"), row.names = c(NA, -270L))

data %>% mutate_if(is.character, as.factor) %>% summary

data %>% filter(type == 'hyacinth', source == 'seed')

cross_flowers('hyacinth 201')

cross_flowers('hyacinth 201 201', single_cross = TRUE)

cross_flowers('hyacinth 001 020 201') # hyacinth available from seed

