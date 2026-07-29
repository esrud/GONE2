# Master vs experimental — visual comparison

Side-by-side runs of the **master** (left) and **experimental**
(`experiments/ga-exploration`, right) branches on the simulation
benchmarks. Each pair was produced by `gone-ncurses` with the
`-f <truth>` overlay, so the same chart shows:

- the GA's Ne estimate (cyan),
- the simulation's ground-truth Ne history (dashed red),
- the best SCval and round counter in the top-right corner box.

Lower SCval = better fit. A faithful curve hugs the red dashed line.
Cases are grouped by demographic shape; sample-size variants (`n=100`
vs `n=10`) are shown beneath each scenario where both exist.

## Constant Ne (`FIX_*`)

Flat ground truth, varied at three population sizes (100, 1 000,
10 000). The easiest case — both branches should track the constant
line.

### `FIX_100`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/fix_100_n100_std.png) | ![](../img/fix_100_n100_exp.png) |

### `FIX_100`, n = 10

| master | experimental |
|:---:|:---:|
| ![](../img/fix_100_n10_std.png) | ![](../img/fix_100_n10_exp.png) |

### `FIX_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/fix_1000_n100_std.png) | ![](../img/fix_1000_n100_exp.png) |

### `FIX_10000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/fix_10000_n100_std.png) | ![](../img/fix_10000_n100_exp.png) |

### `FIX_10000`, n = 10

| master | experimental |
|:---:|:---:|
| ![](../img/fix_10000_n10_std.png) | ![](../img/fix_10000_n10_exp.png) |

## Monotonic change

Smooth Ne trajectory — `INCREASE_*` ramps up over generations,
`DECREASE_*` ramps down. `decreaseb_1000` is a second decrease replicate.

### `INCREASE_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/increase_1000_n100_std.png) | ![](../img/increase_1000_n100_exp.png) |

### `DECREASE_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/decrease_1000_n100_std.png) | ![](../img/decrease_1000_n100_exp.png) |

### `DECREASE_1000`, n = 10

| master | experimental |
|:---:|:---:|
| ![](../img/decrease_1000_n10_std.png) | ![](../img/decrease_1000_n10_exp.png) |

### `DECREASE_1000b`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/decreaseb_1000_n100_std.png) | ![](../img/decreaseb_1000_n100_exp.png) |

## Sudden drop (`DROP_*`)

Ne crashes to a low value and stays there — tests the GA's ability to
follow a sharp downward step.

### `DROP_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/drop_1000_n100_std.png) | ![](../img/drop_1000_n100_exp.png) |

## Recent bottlenecks (`BOT2/3/4`)

A recent dip in Ne followed by recovery. Becomes harder as the
bottleneck moves further back in time (BOT2 most recent, BOT4 oldest)
and as sample size shrinks.

### `BOT2_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/bot2_1000_n100_std.png) | ![](../img/bot2_1000_n100_exp.png) |

### `BOT2_1000`, n = 10

| master | experimental |
|:---:|:---:|
| ![](../img/bot2_1000_n10_std.png) | ![](../img/bot2_1000_n10_exp.png) |

### `BOT3_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/bot3_1000_n100_std.png) | ![](../img/bot3_1000_n100_exp.png) |

### `BOT4_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/bot4_1000_n100_std.png) | ![](../img/bot4_1000_n100_exp.png) |

### `BOT4_1000`, n = 10

| master | experimental |
|:---:|:---:|
| ![](../img/bot4_1000_n10_std.png) | ![](../img/bot4_1000_n10_exp.png) |

## Old bottleneck (`BOTOLD`)

Bottleneck pushed further into the past — the LD signal is weaker.

### `BOTOLD_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/botold_1000_n100_std.png) | ![](../img/botold_1000_n100_exp.png) |

## Soft scenarios (`botsoft`, `expsoft`)

Smoother demographic trajectories — bowl-shaped or expansion.

### `botsoft_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/botsoft_1000_n100_std.png) | ![](../img/botsoft_1000_n100_exp.png) |

### `expsoft_1000`, n = 100

| master | experimental |
|:---:|:---:|
| ![](../img/expsoft_1000_n100_std.png) | ![](../img/expsoft_1000_n100_exp.png) |
