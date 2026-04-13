# data_review_project
这是结课作业,新建的一个公共仓库

# 1,小组成员

1 姓名:欧阳环   学号:2025303120033 github@oyang666

2 姓名:蔡洋洋   学号:2025303120047 github@CYY386

3 姓名:罗秋艳   学号:2025303110041 github@luoqiuqiu123

4 姓名:涂君毅   学号:2025303120036 github@tujunyi124

# 2,研究项目

利用R对文献中的图片进行复现.
这是文献原文:
[2025_全球固氮生物_气候变化-导致多样性下降_Li_NatCommun.pdf](https://github.com/user-attachments/files/26047138/2025_._.-._Li_NatCommun.pdf)

# 3,数据代码来源

翻到文章末尾，在“参考文献 (References)”前面寻找关键信息:

**Data availability (数据可用性)**

**Code availability (代码可用性)**

文献数据源地址:https://doi.org/10.24433/CO.9682287.v2 (此文献数据源提供了云数据复现的功能,能在云端一键运行复现服务)

# 4,运行环境及运行方法

使用R进行本地复现,按照数据处理的源代码,安装必要的包,配置好相关环境,即可开始进行复现.

# 5,复现步骤

1,从https://doi.org/10.24433/CO.9682287.v2 上下载必要的数据.

2,打开Rstudio,安装数据分析所需的包.

3,本地化处理,将源代码中的地址关键字改成本地路径.

4,优化代码,使其更加贴近原图.

5,文后链接中路径是我的已经修改过的路径,如果需要再次复现,请进行修改

# 6,结论

本次进行了Figure1a和Figure1c的绘制,这是作者的原图:

[Figure 1.pdf](https://github.com/user-attachments/files/26047609/Figure.1.pdf)

然后这是利用R复现的图:

<img width="6000" height="3161" alt="Image" src="https://github.com/user-attachments/assets/98b2804c-4271-4f3f-8e3f-6a00874c2b95" />
<img width="6600" height="3480" alt="Image" src="https://github.com/user-attachments/assets/15a55695-7221-43d9-98cb-4a3b234079f7" />

可以看到,在图的主体部分与原图一致,但还是在细节上有点差异,但是作为数据复现来说已经符合要求了.
该研究使用公开宏基因组数据和开源代码进行分析，环境配置简单，因此具有较高的可重复性，适合作为可重复研究案例进行复现分析。

# 7,本次分析所使用代码

这部分是Figure1a的带代码:

```R
# figure1a 
options(warn=-1)
require(data.table)
require(UpSetR)
nf_fixed <- fread("/home/oyang/文档/文献文件/数据复现/sourcedata/Figure S1a.csv")
nf_fixed$'count' <- rep(1,nrow(nf_fixed))
dat2_1 <- dcast(nf_fixed,GCF~sign,fun.aggregate = sum,value.var = c("count"))

test <- unique(nf_fixed[,c(4,7:12,27:28)],by = c("GCF"))
dat2_1 <- dat2_1|>as.data.frame()
dat2_1_0 <- dat2_1
dat2_1_0[,c(2:8)][dat2_1_0[,c(2:8)]>0] <- 1
dat2_1_0 <- test[dat2_1_0,on = c("GCF")]
dat2_1_0$main[dat2_1_0$main == "HDK" & dat2_1_0$`anfG or vnfG`>0] <- "HDGK"
dat2_1_meta <- colSums(dat2_1[,2:8])|>data.frame()
colnames(dat2_1_meta)[1] <- "Genes" 
dat2_1_meta$'sets' <- rownames(dat2_1_meta)
dat2_1_meta <- dat2_1_meta[,c(2,1)]

#require(UpSetR)
p2_1 <- upset(dat2_1_0,
              sets = rev(c("nifH","nifD","nifK","nifE","nifN","nifB","anfG or vnfG")),
              #nsets = 7,
              order.by = c("degree", "freq"),
              decreasing = c(TRUE, TRUE),
              keep.order = TRUE,
              sets.x.label = "Genomes with any nif gene",
              query.legend = "top", 
              point.size = 3,
              line.size = 1,
              mb.ratio = c(0.7,0.3),
              mainbar.y.label = "Genomes",
              sets.bar.color = rev(c("#ffa726",'#42a5f5','#42a5f5',
                                     'grey50','grey50','grey50',"#26c6d5")),
              queries = list(
                list(
                  query = elements,
                  params = list("main","HDK"),
                  color = "#42a5f5",
                  active = T,
                  query.name = "Completed Nitrogenase (nifHDK)"
                ),
                list(
                  query = elements,
                  params = list("main","H"),
                  color = "#ffa726",
                  active = T,
                  query.name = "Incompleted Nitrogenase (nifH only)"
                ),
                list(
                  query = elements,
                  params = list("main","HDGK"),
                  color = "#26c6d5",
                  active = T,
                  query.name = "Incompleted Nitrogenase (anf or vnf)"
                )
              ),
              set.metadata = list(
                data = dat2_1_meta, 
                plots = list(
                  list(
                    type = "hist", 
                    column = "Genes", 
                    colors = rev(c("#ffa726",'#42a5f5','#42a5f5',
                                   'grey50','grey50','grey50',"#26c6d5")),
                    assign = 20
                  ), 
                  list(
                    type = "matrix_rows", 
                    column = "sets", 
                    colors = c(
                      nifH = "#ffa726", 
                      nifD = '#42a5f5',
                      nifK = '#42a5f5',
                      nifE = 'grey50',
                      nifN = 'grey50',
                      nifB = 'grey50',
                      `anfG or vnfG` = "#26c6d5"
                    ),
                    alpha = 0.5
                  )
                )
              )
              
              
)
png("../results/Figure_1a.png", width = 10, height = 5.27, units = "in", res = 600)
print(p2_1) 
dev.off()
```

这部分是Figure1c的代码:

```R
# 加载必要的包
library(data.table)
library(ggplot2)
library(patchwork)
library(readxl)

options(warn = -1)

# 读取数据
file_path <- "/home/oyang/文档/文献文件/数据复现/sourcedata/Figure 1.xlsx"
dt_raw <- as.data.table(read_excel(file_path, sheet = "Figure1c"))

# 查看数据结构和 measure 有哪些值
print(names(dt_raw))
print(unique(dt_raw$measure))

# ---------- 构建 test0（用于总标签：nif_count 和 nif_percent）----------
# 假设 measure 中包含 "nif_count" 和 "nif_percent"
test0 <- dt_raw[measure %in% c("nif_count", "nif_percent")]
# 添加 type 和 label 列
test0[, type := ifelse(measure == "nif_count", "abs", "per")]
test0[, label := ifelse(type == "abs", 
                        as.character(value), 
                        paste0(as.character(round(value, 4) * 100), "%"))]

# ---------- 构建 test1（用于堆叠条形图：no_evid_count, evid_count, no_evid_percent, evid_percent）----------
test1 <- dt_raw[measure %in% c("no_evid_count", "evid_count", "no_evid_percent", "evid_percent")]
test1[, type := ifelse(measure %in% c("no_evid_count", "evid_count"), "abs", "per")]
test1[, level := ifelse(measure %in% c("no_evid_count", "no_evid_percent"), "no_evid", "evid")]
test1[, level := factor(level, levels = c("no_evid", "evid"))]
phylum_levels <- c("Pseudomonadota", "Bacillota",
                   "Thermodesulfobacteriota", "Cyanobacteriota",
                   "Euryarchaeota", "Bacteroidota",
                   "Campylobacterota", "Chlorobiota",
                   "Actinomycetota", "Verrucomicrobiota",
                   "Spirochaetota", "Nitrospirota",
                   "Chloroflexota", "Deferribacterota",
                   "Planctomycetota", "Chrysiogenota",
                   "Kiritimatiellota", "Candidatus Thermoplasmatota",
                   "Aquificota", "Acidobacteriota", "Fusobacteriota")
# 反转顺序（因为 coord_flip 后希望从上到下按列表顺序）
test0[, phylum := factor(phylum, levels = rev(phylum_levels))]
test1[, phylum := factor(phylum, levels = rev(phylum_levels))]

# ---------- 绘图 ----------
p3_1 <- ggplot(test1[type == "abs", ], aes(x = phylum, y = value, fill = level)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = test0[type == "abs", ], aes(label = label, fill = NULL),
            hjust = -0.5, size = 3, color = "black") +
  coord_flip() +facet_grid(. ~ type, 
                           labeller = labeller(type = c(abs = "Nitrogen-fixing genomes", 
                                                        per = "Percentage"))) +
  scale_fill_manual(name = "yes",values = c("#74b1df", "#7594a9")) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        strip.background.x = element_rect(color = "white", fill = "white", size = 1.5, linetype = "solid"),
        strip.placement = "outside",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic", hjust = 0.55),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.05),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 0))

p3_2 <- ggplot(test1[type != "abs", ], aes(y = phylum, x = value, fill = level)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = test0[type != "abs", ], aes(label = label, fill = NULL),
            hjust = 1.3, size = 3, color = "black") +
  scale_fill_manual(name = "No",values = c("#74b1df", "#7594a9")) +
  scale_y_discrete(position = "right") +
  scale_x_reverse(labels = abs, expand = expansion(mult = c(0.1, 0))) +
  facet_grid(. ~ type, 
             labeller = labeller(type = c(abs = "Nitrogen-fixing genomes", 
                                          per = "Percentage"))) +
  theme_classic() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        strip.background.x = element_rect(color = "white", fill = "white", size = 1.5, linetype = "solid"),
        strip.placement = "outside",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 0, b = 10, l = 10))

p <- p3_2 | p3_1

# 保存图片
png("/home/oyang/文档/文献文件/数据复现/review/Figure_1c.png", width = 11, height = 5.8, units = "in", res = 600)
print(p)
dev.off()
```

