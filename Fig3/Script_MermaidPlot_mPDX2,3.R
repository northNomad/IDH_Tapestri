wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(DiagrammeR)
library(scales)
library(RColorBrewer)
library(DiagrammeRsvg)
library(rsvg)
library(ggsci)

library(scales)
pal_simpsons()(6)
show_col(pal_simpsons()(6))

r <- c(94, 735, 2024, 867, 3156)
r <- r / sum(r)
r*10

#F2CB9E, #FEE7D2, #F1F5F8, #D5E4F4, #9DCBE4, #A8C8EB
grViz("
digraph a_nice_graph {


# node [fontname=Helvetica, shape=ellipse, style=filled]
## Define roots
node [shape = plaintext, fontname=Helvetica]


PM160345_root [label='root-PM160345', color='#FFF7F3']
PM150437_root [label='root-PM150437', color='#FFF7F3']




node [shape=circle, stroke=black, fontname=Helvetica, style=filled, fixedsize=TRUE]

### PM160345
PM160345_DNMT3A [label='<I>DNMT3A</I>@^{R882C}<br/>N=94', color='#FED439FF', width=0.1335986]
PM160345_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=735', color='#FED439FF', width=1.0446276]
PM160345_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=2024', color='#FED439FF', width=2.8766345]
PM160345_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=867', color='#FED439FF', width=1.2322342]

PM160345_root -> PM160345_DNMT3A
PM160345_DNMT3A -> PM160345_IDH1
PM160345_IDH1 -> PM160345_NPM1
PM160345_NPM1 -> PM160345_FLT3

### PM150437
PM150437_TP53 [label='<I>TP53</I>@^{I255S}<br/>N=3156', color=DarkTurquoise, width=4.5898778]

PM150437_root -> PM150437_TP53

rankdir=LR
}

[1]: 'left'
[2]: 10:20
") -> p_scale_across_pdxc
p_scale_across_pdxc

p_scale_across_pdxc %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg_svg("PDXC_p_scale_across_samples.svg")

# p_scale_by_sample %>%
#   DiagrammeRsvg::export_svg() %>%
#   charToRaw() %>%
#   rsvg_svg("PDXA_p_scale_by_sample.svg")

grViz("
digraph a_nice_graph {


# node [fontname=Helvetica, shape=ellipse, style=filled]
## Define roots
node [shape = plaintext, fontname=Helvetica]


PM160345_root [label='root-PM160345', color='#FFF7F3']
PM150437_root [label='root-PM150437', color='#FFF7F3']




node [shape=ellipse, stroke=black, fontname=Helvetica, style=filled]

### PM160345
PM160345_DNMT3A [label='<I>DNMT3A</I>@^{R882C}<br/>N=94', color='#FED439FF']
PM160345_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=735', color='#FED439FF']
PM160345_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=2024', color='#FED439FF']
PM160345_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=867', color='#FED439FF']

PM160345_root -> PM160345_DNMT3A
PM160345_DNMT3A -> PM160345_IDH1
PM160345_IDH1 -> PM160345_NPM1
PM160345_NPM1 -> PM160345_FLT3

### PM150437
PM150437_TP53 [label='<I>TP53</I>@^{I255S}<br/>N=3156', color=DarkTurquoise]

PM150437_root -> PM150437_TP53

rankdir=LR
}

[1]: 'left'
[2]: 10:20
") -> p_no_scale_pdxc
p_no_scale_pdxc

p_no_scale_pdxc %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg_svg("PDXC_p_no_scale.svg")


# PDXD --------------------------------------------------------------------

r_d <- c(30, 409, 2035, 572, 151, 543, 124)
r_d <- r_d / sum(r_d)
r_d * 10


grViz("
digraph a_nice_graph {


# node [fontname=Helvetica, shape=ellipse, style=filled]
## Define roots
node [shape = plaintext, fontname=Helvetica]


PM160345_root [label='root-PM160345', color='#FFF7F3']
PM160053_root [label='root-PM160053', color='#FFF7F3']


node [shape=circle, stroke=black, fontname=Helvetica, style=filled, fixedsize=TRUE]

### PM160345
PM160345_DNMT3A [label='<I>DNMT3A</I>@^{R882C}<br/>N=30', color='#FED439FF', width=0.07763975]
PM160345_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=409', color='#FED439FF', width=1.05848861]
PM160345_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=2035', color='#FED439FF', width=5.26656315]
PM160345_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=572', color='#FED439FF', width=1.48033126]

PM160345_root -> PM160345_DNMT3A
PM160345_DNMT3A -> PM160345_IDH1
PM160345_IDH1 -> PM160345_NPM1
PM160345_NPM1 -> PM160345_FLT3

### PM160053
PM160053_IDH2 [label='<I>IDH2</I>@^{R140Q}<br/>N=151', color=DarkTurquoise, width=0.39078675]
PM160053_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=543', color=DarkTurquoise, width=1.40527950]
PM160053_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=124', color=DarkTurquoise, width=0.32091097]

PM160053_root -> PM160053_IDH2
PM160053_IDH2 -> PM160053_NPM1
PM160053_NPM1 -> PM160053_FLT3

rankdir=LR
}

[1]: 'left'
[2]: 10:20
") -> p_scale_across_pdxd
p_scale_across_pdxd

p_scale_across_pdxd %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg_svg("PDXD_p_scale_across_samples.svg")
