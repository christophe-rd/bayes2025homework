# Homework 1: 
## Assigned 15 January 2025

## This homework has four parts:

1 -- Simulate data for a linear regression problem (continuous predictor and response variables) of your choice and show a linear model can return the parameters you set. Your n should be 100 and you should set your error (sigma) to a level you would expect in an ecological or evolutionary biology study. 

Next, add a binary covariate (for example, a treatment and control) to your dataset that interacts with your other predictor variable. Make the effect size of this interaction HALF the effect size of the first predictor variable, and again test how well you can return the parameter values.

Step 2 -- Building on the above example, let's look at how the estimates of parameters returned from your linear model change with sample size, and make our analyses a little more robust. 

Using just 20% of your data, see how well you can return the parameters you set.

Next, let's make a plot of sampling from your data from a low sample size (e.g., 2 or 5) to the n of your data (100) showing the estimated parameters for each sample size (you do not need to sample by 1, so you can do, 5, 15 etc.). Make a plot with your 1:100 on the horizontal axis and your estimated parameters on the vertical (you need either as many plots as parameters or a way to show all the different parameters
 on one plot). Compare how well the model does across the different parameters. Which is better or worse at and why? 

[comment]: <> (Some people struggled on what I meant here, I updated it after the assignment but could use more work.)
Let's improve our sampling now. So far we have take just ONE draw from our set of parameters which means our sigma term has some Monte Carlo error in it, so let's take 10 draws each time (so we need to set up a loop or such in R that samples from 1 (or 2 or 5 etc.) to 100 and _at each step it repeats simulating the data and getting the estimates 10 times_). Re-make your plot. (If you get stuck here for a while, don't panic, but move onto Steps 3-4.)

Step 3 -- Now, repeat the above but change your error term by 50%. Compare these results to what you found before and explain why you think they changed. 

**Challenge Options** (optionals extra tasks to add if you finish this quickly and are unsurprised by the results): 

Challenge A -- Make your treatment have three levels (control, treatment 1, treatment 2) and use the dummy variable approach to repeat the above steps. 

Challenge B -- Add an additional binary treatment (so you have treatment 1 which has control and treatment levels and you add a treatment 2 which also has control and treatment levels) in a full factorial design with treatment 1 and repeat the above steps including all two-way and three-way interactions. 

Challenge C -- Repeat step 3, but also extract an estimate of uncertainty from your posterior (e.g., 50% uncertainty intervals) and plot the results. 

4 -- Review the [sample datasets](https://github.com/temporalecologylab/bayes2025homework/blob/main/analyses/datasetsbayes.R) and pick one you want to use to go through the complete workflow with (which is the homework for next week). Ideally, you'll get to work on that one, but I may ask people to pick a new one so everyone is not using the same one. 


## How long should my answers be?
Your answers should be BRIEF. We'll review the homework in class on Tuesday. 

### What file types can I submit?
I prefer .R or .Rmd for code, but Quarto works too. If you use R then you can write answers with comments or submit a separate .md or .txt file with your answers. Please include all requested plots and complete code (do overwrite code to do the last step). Compile your qmd and md files and submit the html or pdf with the .md and .qmd.

(Do you happen to be writing your code in Word or Google Docs? Don't worry, many good scientists have done this, but now would be the time to stop.)

### What models should I fit?
Use rstanarm (`stan_glm`). Either way, use the mean point estimates for this exercise. If your model does not converge **do not report the results** just record `NA` for that output. Rstanarm vignette [here](https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html).

# Due 20 January 2025 by 3pm

### Submit how? If you are enrolled in class, you must submit homework.
You can submit via Canvas or GitHub. To submit on GitHub:

1. Make a folder (name it your first name, all lowercase letters) in homeworksubmitted/

2. Put your analyses and any output there, format the code so it runs from the repo without additional files.

3. Submit this as [a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request). 

If you submit on Canvas, follow steps 1-2, then compress your folder and upload it via the assignment. 

If you're auditing you don't have to submit your homework, but I still highly recommend doing it (it will test your skills in the course and make it easier for you to enjoy the homework review in class each week). You're welcome to submit your homework to me if you are auditing the course though please try to submit via one of the two methods above (I can add you to Canvas if you are auditing but need enough notice), because my email inbox is not a great place for homework.

