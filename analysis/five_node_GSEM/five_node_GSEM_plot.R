library(dplyr)
library(tidyr)
library(ggplot2)


merged_results <- read.csv('/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_GSEM_eval/nesmr_vs_genomic_sem_results.csv') %>%
    mutate(bias = beta - true_beta) %>%
    separate(edge, c('from', 'to'), sep = '->')

merged_results %>%
    group_by(from, to, model, n) %>%
    summarise(
        bias = mean(bias),
        max_bias = max(abs(bias))
    ) %>%
    arrange(desc(abs(max_bias)))

p1 <- filter(merged_results, n == 40000) %>% ggplot(.) +
    facet_grid(rows = vars(from), cols = vars(to)) +
    xlim(c(-0.5, 0.5)) +
    geom_histogram(aes(x = bias, fill = model), alpha = 0.5, position="identity") +
    ggtitle('N = 40,000') +
    ylab('') +
    theme_minimal(base_size = 24) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- filter(merged_results, n == 10000) %>% ggplot(.) +
    facet_grid(rows = vars(from), cols = vars(to)) +
    xlim(c(-0.5, 0.5)) +
    geom_histogram(aes(x = bias, fill = model), alpha = 0.5, position="identity") +
    ggtitle('N = 10,000') +
     ylab('') +
    theme_minimal(base_size = 24) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave('nesmr_vs_genomic_sem_bias_10000.jpeg', p2, width = 14, height = 14, dpi = 300)
ggsave('nesmr_vs_genomic_sem_bias_40000.jpeg', p1, width = 14, height = 14, dpi = 300)
