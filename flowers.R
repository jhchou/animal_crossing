library(tidyverse)

data <- read.csv(file = 'flowers.csv', stringsAsFactors = FALSE) %>% # type + color will be factors
  as_tibble() %>%
  rename_all(tolower) %>%
  mutate(
    source = as.factor(ifelse(source == 1, 'Seed', ifelse(source == 2, 'Island', '(cross)'))),
    seedbag = (seedbag != 0),
    breedtrue = ((red != 1) & (yellow != 1) & (white != 1) & (brightness != 1))
  ) %>%
  select(-id, -seedbag)

data %>%
  select(type, color) %>%
  unique %>%
  group_by(type) %>%
  summarise(num_color = n()) %>%
  arrange(desc(num_color))

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

cross_flower <- function(type, r1, y1, w1, r2, y2, w2, b1 = 0, b2 = 0) {
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
    group_by(type,red,yellow,white,brightness, color, breedtrue) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    group_by(color) %>%
    mutate(prop = n / sum(n)) %>%
    select(-n) %>%
    arrange(color, desc(prop))
}

cross_flower('Hyacinth',1,0,1,1,0,1)

