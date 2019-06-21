# Tissue-specificity metrics

tspex provides twelve distinct tissue-specificity metrics, which differ in their assumptions, scale and properties. Broadly, these metrics can be divided into two groups[^1]:

- **General scoring metrics:** Summarize in a single value how tissue-specific or ubiquitous is a gene across all tissues.
- **Individualized scoring metrics:** Quantify how specific is the expression of each gene to each tissue.

!!! warning ""
    For the following equations $x_i$ represents the gene expression in tissue $i$ and $n$ is the number of tissues. Some metrics give results that are not in the $[0,1]$ interval, such as the Simpson index, that lies in the $[\frac{1}{n},1]$ range. For these metrics, we provide transformed versions, denoted by the $'$ symbol, that range from 0 (perfectly ubiquitous) to 1 (perfectly tissue-specific).

## General scoring metrics

The general scoring metrics included in tspex are: Counts[^2], Tau[^3], Gini coefficient[^4], Simpson index[^5], Shannon entropy specificity[^6], ROKU specificity[^7], Specificity measure dispersion (SPM DPM)[^8], and Jensen-Shannon specificity dispersion (JSS DPM).

### Counts

$$
\begin{align}
    & t(x_i) =
    \begin{cases}
        1, & \text{if } x_i > threshold \\
        0, & \text{if } x_i \le threshold
    \end{cases} \\[9pt]
    & \operatorname{Counts} = \frac{n-\sum_{i=1}^n(t(x_i))}{n-1}
\end{align}
$$

### Tau

$$
\begin{align}
    & \widehat{x_i} = \frac{x_i}{\max\limits_{0\le i \le n}(x_i)} \\[9pt]
    & \operatorname{Tau} = \frac{\sum_{i=1}^n(1-\widehat{x_i})}{(n-1)}
\end{align}
$$

### Gini coefficient

$$
\begin{align}
    & X = (x_1,x_2,\ldots,x_i,\ldots,x_{n-1},x_{n}) \\[9pt]
    & \qquad \;\; x_1 \le x_2 \le \ldots x_i \le \ldots \le x_{n-1} \le x_{n} \notag\\[9pt]
    & \operatorname{Gini} = \frac{\sum_{i=1}^n(2i-n-1)x_i}{n\sum_{i=1}^n(x_i)} \\[9pt]
    & \operatorname{Gini'} = \operatorname{Gini} \frac{n}{n-1}
\end{align}
$$

### Simpson index

$$
\begin{align}
    & p_i = \frac{x_i}{\sum_{i=1}^n(x_i)} \\[9pt]
    & \operatorname{Simpson} = \sum_{i=1}^n(p_i^2) \\[9pt]
    & \operatorname{Simpson'} = \frac{\operatorname{Simpson} - \frac{1}{n}}{1 - \frac{1}{n}}
\end{align}
$$

### Shannon entropy specificity (HS)

$$
\begin{align}
    & p_i = \frac{x_i}{\sum_{i=1}^n(x_i)} \\[9pt]
    & H = -\sum_{i=1}^n(p_i \log_2(p_i)) \\[9pt]
    & \operatorname{HS} = log_2(n) - H \\[9pt]
    & \operatorname{HS'} = \frac{\operatorname{HS}}{log_2(n)}
\end{align}
$$

### ROKU specificity

$$
\begin{align}
    & M = \operatorname{median}_{0\le i \le n}(x_i) \\[9pt]
    & S = \operatorname{median}_{0\le i \le n}(|x_i - M|) \\[9pt]
    & u_i = \frac{x_i - M}{5S + 10^{-4}} \\[9pt]
    & w(u_i) =
    \begin{cases}
        (1-u^2)^2, & \text{if } |u_i| \le 1 \\
        0, & \text{if } |u_i| > 1
    \end{cases} \\[9pt]
    & t = \frac{\sum_{i=1}^n(x_i w(u_i))}{\sum_{i=1}^n(w(u_i))} \\[9pt]
    & X' = (|x_1-t|,|x_2-t|,\ldots,|x_i-t|,\ldots,|x_{n-1}-t|,|x_{n}-t|) \\[9pt]
    & p_i = \frac{x_i'}{\sum_{i=1}^n(x_i')} \\[9pt]
    & H = -\sum_{i=1}^n(p_i \log_2(p_i)) \\[9pt]
    & \operatorname{ROKU} = log_2(n) - H \\[9pt]
    & \operatorname{ROKU'} = \frac{\operatorname{ROKU}}{log_2(n)}
\end{align}
$$

### Specificity measure dispersion (SPM DPM)

$$
\begin{align}
    & \operatorname{\overline{SPM}} = \frac{\sum_{i=1}^n(\operatorname{SPM_i})}{n} \\[9pt]
    & \sigma = \sqrt{\frac{\sum_{i=1}^n(\operatorname{SPM_i} - \operatorname{\overline{SPM}})^2}{n-1}} \\[9pt]
    & \operatorname{SPM\ DPM} = \sigma \sqrt{n}
\end{align}
$$

### Jensen-Shannon specificity dispersion (JSS DPM)

$$
\begin{align}
    & \operatorname{\overline{JSS}} = \frac{\sum_{i=1}^n(\operatorname{JSS_i})}{n} \\[9pt]
    & \sigma = \sqrt{\frac{\sum_{i=1}^n(\operatorname{JSS_i} - \operatorname{\overline{JSS}})^2}{n-1}} \\[9pt]
    & \operatorname{JSS\ DPM} = \sigma \sqrt{n}
\end{align}
$$

## Individualized scoring metrics

The general scoring metrics included in tspex are: Tissue-specificity index (TSI)[^9], Z-score[^10], Specificity measure (SPM)[^11], and Jensen-Shannon specificity (JSS)[^12].

### Tissue-specificity index (TSI)

$$
\begin{align}
    & \operatorname{TSI_i} = \frac{x_i}{\sum_{i=1}^n(x_i)}
\end{align}
$$

### Z-score

$$
\begin{align}
    & \overline{x} = \frac{\sum_{i=1}^n(x_i)}{n} \\[9pt]
    & \sigma = \sqrt{\frac{\sum_{i=1}^n(x_i - \overline{x})^2}{n-1}} \\[9pt]
    & \operatorname{Z-Score_i} = \frac{x_i - \overline{x}}{\sigma} \\[9pt]
    & \operatorname{Z-Score_i'} = \frac{\operatorname{Z-Score_i} + \frac{(n-1)}{\sqrt{n}}}{2 \frac{(n-1)}{\sqrt{n}}}
\end{align}
$$

### Specificity measure (SPM)

$$
\begin{align}
    & X = (x_1,x_2,\ldots,x_i,\ldots,x_{n-1},x_{n}) \\[9pt]
    & \operatorname{SPM_i} = \frac{x_i^2}{x_i*||X||_2}
\end{align}
$$

### Jensen-Shannon specificity (JSS)

$$
\begin{align}
    & X = (x_1,x_2,\ldots,x_i,\ldots,x_{n-1},x_{n}) \\[9pt]
    & X' = (0,0,\ldots,x_i,\ldots,0,0) \\[9pt]
    & p_i = \frac{x_i}{\sum_{i=1}^n(x_i)} \\[9pt]
    & q_i = \frac{x'_i}{\sum_{i=1}^n(x'_i)} \\[9pt]
    & \operatorname{JSS_i} = 1 - \sqrt{\frac{\sum_{i=1}^n(p_i \log_2(p_i))}{2} - \sum_{i=1}^n\left(\frac{p_i+q_i}{2} \log_2\left(\frac{p_i+q_i}{2}\right)\right)}
\end{align}
$$

## References

[^1]: Kryuchkova-Mostacci, Nadezda, and Marc Robinson-Rechavi. "A benchmark of gene expression tissue-specificity metrics." *Briefings in bioinformatics* 18.2 (2017): 205-214.
[^2]: Duret, Laurent, and Dominique Mouchiroud. "Determinants of substitution rates in mammalian genes: expression pattern affects selection intensity but not mutation rate." *Molecular biology and evolution* 17.1 (2000): 68-070.
[^3]: Yanai, Itai, et al. "Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification." *Bioinformatics* 21.5 (2004): 650-659.
[^4]: Ceriani, Lidia, and Paolo Verme. "The origins of the Gini index: extracts from Variabilità e Mutabilità (1912) by Corrado Gini." *The Journal of Economic Inequality* 10.3 (2012): 421-443.
[^5]: Simpson, Edward H. "Measurement of diversity." *Nature* 163.4148 (1949): 688.
[^6]: Schug, Jonathan, et al. "Promoter features related to tissue specificity as measured by Shannon entropy." *Genome biology* 6.4 (2005): R33.
[^7]: Kadota, Koji, et al. "ROKU: a novel method for identification of tissue-specific genes." *BMC bioinformatics* 7.1 (2006): 294.
[^8]: Pan, Jian-Bo, et al. "PaGeFinder: quantitative identification of spatiotemporal pattern genes." *Bioinformatics* 28.11 (2012): 1544-1545.
[^9]: Julien, Philippe, et al. "Mechanisms and evolutionary patterns of mammalian and avian dosage compensation." *PLoS biology* 10.5 (2012): e1001328.
[^10]: Vandenbon, Alexis, and Kenta Nakai. "Modeling tissue-specific structural patterns in human and mouse promoters." *Nucleic acids research* 38.1 (2009): 17-25.
[^11]: Xiao, Sheng-Jian, et al. "TiSGeD: a database for tissue-specific genes." *Bioinformatics* 26.9 (2010): 1273-1275.
[^12]: Cabili, Moran N., et al. "Integrative annotation of human large intergenic noncoding RNAs reveals global properties and specific subclasses." *Genes & development* 25.18 (2011): 1915-1927.