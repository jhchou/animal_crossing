library(tidyverse)

data <- read.csv(file = 'flowers.csv', stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename_all(tolower) %>%
  mutate(
    source = ifelse(source == 1, 'seed', ifelse(source == 2, 'island', 'cross')),
    seedbag = (seedbag != 0),
    breedtrue = ((red != 1) & (yellow != 1) & (white != 1) & (brightness != 1))
  ) %>%
  mutate_if(is.character, tolower) %>%
  select(-id, -seedbag)

cross_flower_full <- function(type, r1, y1, w1, b1, r2, y2, w2, b2) {
  # for all except rose, b1 and b2 = 0

  cross_genotypes <- function(g1, g2) {
    # g1 / g2 genotypes: either 0, 1, or 2
    # return possible genotypes
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

  expand.grid(
    type = type,
    red = cross_genotypes(r1, r2),
    yellow = cross_genotypes(y1, y2),
    white = cross_genotypes(w1, w2),
    brightness = cross_genotypes(b1, b2),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    left_join(data %>% select(-source), by = c("type", "red", "yellow", "white", "brightness")) %>%
    mutate(total_n = n()) %>%
    group_by(type, red, yellow, white, brightness, color, breedtrue, total_n) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    group_by(color) %>%
    mutate(prop_cross = n/ total_n, prop_color = n / sum(n)) %>%
    ungroup() %>%
    select(-n, -total_n) %>%
    arrange(color, desc(prop_color))
}
# cross_flower_full('hyacinth',1,0,1,0,1,0,1,0)

cross_flower <- function(input) {
  temp <- unlist(strsplit(input, ' '))
  type <- tolower(temp[1])
  genotype1 <- temp[2]
  genotype2 <- temp[3]
  
  g1 <- as.numeric(unlist(strsplit(genotype1, split = '')))
  g2 <- as.numeric(unlist(strsplit(genotype2, split = '')))
  r1 = g1[1]
  y1 = g1[2]
  w1 = g1[3]
  b1 = ifelse(length(g1) == 4, g1[4], 0)
  r2 = g2[1]
  y2 = g2[2]
  w2 = g2[3]
  b2 = ifelse(length(g2) == 4, g2[4], 0)
  cross_flower_full(type, r1, y1, w1, b1, r2, y2, w2, b2)
}

cross_flower('hyacinth 101 101')
cross_flower('rose 2020 1010')
