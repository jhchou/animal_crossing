library(tidyverse)

# Read flower genotype / phenotype / source data
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
# - cross_flower('tulip 001 020')
# - b1 and b2 = 0 for all except rose
cross_flower <- function(input) {
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


cross_flower('hyacinth 101 101')
cross_flower('rose 2020 1010')


cross_flower_all <- function(type, genotypes) {
  # Input any number of genotypes and generate results of all possible crosses
  # Combinations with repetitions: https://stackoverflow.com/questions/55738847/r-creating-combinations-with-replacement
  type <- tolower(type)
  genotypes <- unlist(strsplit(genotypes, ' '))
  
  result <- expand.grid(genotypes, genotypes, stringsAsFactors = FALSE) %>%
    mutate(key = paste(pmin(Var1, Var2), pmax(Var1, Var2), sep = "-")) %>%
    filter(!duplicated(key)) %>%
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
    mutate(output = map(input, cross_flower)) %>%
    unnest(output) %>%
    select(cross, genotype, color, breedtrue, prop_cross, prop_color)
  
  return(result)
}

cross_flower_all('hyacinth', '001 222')
