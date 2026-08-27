---
layout: post
title: "Propensity Scores: What They Do, and What They Don't"
date: 2026-08-27
tags: [causal inference, propensity scores]
math: true
---

Start with a story.

Older patients are more likely to receive the treatment; younger patients mostly end up in the control group. But younger patients recover better than older ones no matter what you give them. So even if the treatment genuinely helps, a naive comparison of outcomes between the treated and the untreated can make the treatment look *worse* than doing nothing. Age is doing the work, not the drug.

That is confounding, and it is the problem propensity scores were built to attack. What follows is what a propensity score is, the assumptions that make it work, the single most common way people misuse it, and the two main ways to actually use one.

## The counterfactual setup

For each person we would like to know two things: their outcome if treated, and their outcome if untreated. The average causal effect is the difference between those two means across the population.

The difficulty is structural. Each person contributes exactly one of the two potential outcomes to the data. The other is missing. Everything below is machinery for filling that gap with people who resemble the ones we cannot observe.

Confounders are the covariates associated with both treatment assignment and the outcome. They are why the observed association is not the causal effect — exactly as in the story above, where age drives both who gets treated and who recovers.

## What the propensity score is

Define

$$\pi(L) = \Pr[A = 1 \mid L]$$

the conditional probability of receiving treatment given the measured covariates $L$. This is a person's *propensity* to receive treatment given what you measured about them, which is where the name comes from.

In an ideal randomized trial with 1:1 allocation, $\pi(L) = 0.5$ for everyone, whatever $L$ is. In observational data, treatment assignment is outside the investigator's control, so the true $\pi(L)$ varies across people and is unknown. It has to be estimated — usually by a logistic regression of $A$ on the confounders.

## The three assumptions

An observational study can be conceptualized as a conditionally randomized experiment — one in which the probability of assignment depends on the covariates $L$, rather than being the same for everyone — only if three conditions hold.

**Consistency.** The observed outcome for a treated person equals the outcome she would have had under treatment, and the observed outcome for an untreated person equals the outcome she would have had if she had remained untreated. Implicit in this is *no interference*: whether one person receives treatment has no effect on anyone else's potential outcomes.

**Exchangeability.** The probability of receiving each treatment level depends only on the measured covariates $L$ — no unmeasured confounders. Put plainly: after adjusting for $L$, any remaining outcome difference between treated and untreated is attributable to treatment itself.

**Positivity.** Everyone must have had some chance of receiving each treatment level, i.e. no one has a propensity score of exactly 0 or 1.

Positivity fails in two very different ways:

- **Structural violations.** Certain covariate values make treatment impossible. The canonical example: when estimating the effect of an occupational chemical exposure, being off work is a confounder *and* makes exposure impossible, so the probability of treatment among people off work is exactly zero. Here you cannot recover the effect in the whole population at all — inference has to be restricted to the region where positivity holds.
- **Random violations.** Your sample is finite, so once you stratify on several confounders you will find empty cells by chance. That zero is an artifact of sample size, not structure. When you fit a logistic model for $\Pr[A=1 \mid L]$, the model smooths over those cells by interpolating from everyone else.

The consequence is worth stating plainly: every time you use parametric estimation of the propensity score in the presence of zero cells, you are effectively assuming random nonpositivity.

## The mistake

There is a persistent belief that the treatment model — the model you use to estimate $\Pr[A = 1 \mid L]$ — should predict $A$ as well as possible. It shouldn't.

Take the belief to its limit. If you predicted treatment perfectly, every treated person would have $\pi(L) = 1$ and every untreated person $\pi(L) = 0$. There would be no overlap at all: positivity fails completely.

A prediction mindset also invites *self-inflicted bias*. Colliders are often excellent predictors, and conditioning on them creates systematic bias. Instruments — variables that affect treatment but not the outcome except through treatment — can amplify bias from unmeasured confounders when included.

So: propensity models do not need to predict treatment well. They need to contain the confounders that guarantee exchangeability.

## Using the propensity score, option 1: inverse probability weighting

IP weighting asks a specific question: what if we reweighted people so that treatment no longer depended on the covariates at all — so that assignment became as good as randomized?

Weight each person by the inverse of the conditional probability of receiving the treatment level they actually received. For binary treatment, a treated person with covariates $L$ gets $1 / \Pr[A=1 \mid L]$ and an untreated person gets $1 / \Pr[A=0 \mid L]$.

> A caution that trips people up: the denominator of the IP weight is not $\pi(L)$ for everybody. It is $\pi(L)$ for the treated and $1 - \pi(L)$ for the untreated.

The result is a **pseudo-population** in which $A$ and $L$ are independent. There is no confounding by $L$ there. You fit the plain associational model

$$E[Y \mid A] = \theta_0 + \theta_1 A$$

by weighted least squares in the pseudo-population, and $\hat\theta_1$ consistently estimates the average causal effect. The weighting is what turns an associational regression into a causal one.

**Missing outcomes: the same machinery, applied twice.** People with missing outcomes are censored, and analyzing only the uncensored can induce selection bias. The fix is a second weight. Let $C = 1$ denote censored, and set

$$W^C = \frac{1}{\Pr[C = 0 \mid L, A]} \text{ for the uncensored}, \qquad W^C = 0 \text{ for the censored}$$

then use $W^{A,C} = W^A \times W^C$. The censoring weight upweights the uncensored so they stand in for the similar people who dropped out, reconstructing a pseudo-population in which essentially no one is censored.

**Variance.** Because the weights are themselves estimated, the naive standard error is wrong. Two practical options:

1. **Nonparametric bootstrap.** It propagates uncertainty from *both* stages — estimating the propensity model and fitting the structural model. The cost is computation: it requires real resources, or lots of patience, on large databases.
2. **The robust (sandwich) variance estimator**, standard in most software. Valid but **conservative** — it covers the true parameter more than 95% of the time, so intervals are wider than they need to be.

(A third option, deriving the analytic variance, is correct and efficient but usually means programming the estimator yourself.)

## Using the propensity score, option 2: matching

Pair each treated individual with one or more untreated individuals having a close value of $\pi(L)$ — within 0.05, say. The matched population has, by construction, similar score distributions in both arms.

Choosing how close is close entails a bias–variance tradeoff. Loose criteria match dissimilar people and exchangeability fails; tight criteria discard many people and the intervals widen.

### Check that it worked

Matching is a *procedure*, not a guarantee. Balance can still fail if the model used to estimate the propensity score was misspecified, so the balance diagnostics are not optional — and in practice they are frequently skipped.

The workhorse statistic is the **standardized mean difference** (SMD): the between-group difference in means scaled by the pooled standard deviation,

$$\mathrm{SMD} = \frac{\bar{X}_1 - \bar{X}_2}{\sqrt{(S_1^2 + S_2^2)/2}}$$

with an analogous form for binary covariates using $\hat p(1-\hat p)$ in place of the variance. An absolute SMD above 0.1 is the conventional signal of imbalance. If the SMDs come back above threshold, the response is to respecify the propensity model — add splines for continuous covariates, add interactions — not to accept the imbalance and move on.

### Who are you estimating the effect for?

The deeper problem with matching is not balance but the estimand. Matching each treated person to untreated controls and dropping the unmatched untreated targets the effect *in the treated*. In practice, though, some treated people also fail to find a match, so the estimand becomes a hard-to-describe subset of the population defined by which score values happened to find partners.

That matching automatically enforces positivity is often sold as a strength. There is a price, and it is paid in transportability. Suppose you conclude you can only estimate the effect among people with $\hat\pi(L) < 0.67$. Who are those people? Nobody walks around with a propensity score tattooed on their forehead, and the same score value means different things in different settings. Restricting the study population to the overlapping range is a lazy way to ensure positivity: it buys you a valid estimate for a population you cannot describe.

Restricting on real-world variables instead is more work and much more useful. In the NHEFS smoking-cessation example that runs through *What If*, the two treated individuals with estimated scores above 0.67 turned out to be precisely the people over age 50 who had smoked for fewer than 10 years. Excluding them *by those variables*, and saying so, yields a target population a reader can actually reason about.

## Outcome regression, and having it both ways

The alternative to modeling treatment is modeling the outcome. Fit a model for $E[Y \mid A, C=0, L]$, use it to predict each person's outcome under treatment and under control, and average the difference across the sample. This uses the covariates directly rather than through a score.

The two families lean on different assumptions, and that is the useful part. IP weighting needs the *treatment* model right; outcome regression needs the *outcome* model right. Which means:

- **Run both.** Large disagreement between the IP-weighted and outcome-regression estimates is a warning that at least one model is seriously misspecified. Agreement is reassuring but not proof — both could be wrong in the same direction.
- **Better, use a doubly robust estimator**, which combines a treatment model and an outcome model in a single estimator. Under exchangeability and positivity, it is consistent if *either* model is correct, and you do not need to know which one. Two chances to get it right.

## A working checklist

1. **Justify the confounder set on subject-matter grounds.** Resist the urge to maximize prediction of treatment. Drop instruments. Drop colliders. Keep confounders.
2. **Fit the treatment model** for $\Pr[A=1 \mid L]$ to get a score for each subject. Allow nonlinearities and interactions using splines when justified.
3. **Plot $\hat\pi(L)$ by treatment arm.** Look for overlap. Observed non-overlap is not automatically a structural positivity violation — in a finite sample it may simply be a random one. Decide which you are facing: structural (restrict the population, and describe the restriction in real-world variables) or random (let the model interpolate, and say so).
4. **Then, depending on the method:**
   - *Matching:* check covariate balance in the matched sample, and treat an absolute SMD above 0.1 as a signal to go back to step 2 rather than proceed; report how many treated and untreated subjects were dropped; state the estimand you are left with in substantive terms, not as a range of scores.
   - *Weighting:* assess covariate balance after weighting by comparing the weighted treated and untreated groups using weighted SMDs; deal with extreme weights by truncation or trimming; add censoring weights if outcomes are missing.
5. **Get the uncertainty right.** For weighting, robust or bootstrap variance — the robust estimator is conservative, the bootstrap is expensive. For matching, use a variance estimator that accounts for the matched design.
6. **Estimate the effect a second way** and compare, or use a doubly robust estimator from the start.

## Closing thought

A propensity score analysis lets you conceptualize an observational study as a conditionally randomized experiment and estimate causal effects from it — which matters, because randomized trials are expensive and often difficult to run. What it cannot do is manufacture exchangeability, conjure positivity where none exists, or substitute for knowing your subject matter well enough to say which variables belong in $L$.

---

### References

- Hernán, M. A., & Robins, J. M. (2020). *Causal Inference: What If*. Boca Raton: Chapman & Hall/CRC. Chapters 12, 13, and 15. Available free at [the authors' page](https://miguelhernan.org/whatifbook).
- Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the propensity score in observational studies for causal effects. *Biometrika*, 70(1), 41–55.
- Zhang, Z., Kim, H. J., Lonjon, G., & Zhu, Y. (2019). Balance diagnostics after propensity score matching. *Annals of Translational Medicine*, 7(1), 16.
- Austin, P. C. (2009). Balance diagnostics for comparing the distribution of baseline covariates between treatment groups in propensity-score matched samples. *Statistics in Medicine*, 28(25), 3083–3107.
- Bang, H., & Robins, J. M. (2005). Doubly robust estimation in missing data and causal inference models. *Biometrics*, 61(4), 962–973.
