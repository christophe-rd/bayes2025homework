# Homework 1: 
## Assigned 22 January 2025

## The goal of this homework is to try to complete the workflow for one (given) dataset


As a reminder the steps of the workflow we're using are:

1. Come up with a model.

2. Simulate test data for this model to both make sure you can get your set parameters back from your model (effectively, test your code) and understand better how your model performs. 

3. Set your priors and run some checks on them to make sure they conform with reality and you feel good about them.

4. Run your model on your empirical data. 

5. Perform retrodictive checks using the model fit to your empirical data and consider if that makes you want to change your model.

Remember that you need well-running models at step 2 and 4 (no divergences, Rhat close to 1, good Neffective).

## For this homework, please complete steps 1, 2 and 4 of the workflow on your chosen dataset.
### Please also TRY to complete ONE of the remaining steps (3 or 5).

You must use one of the three [sample datasets](https://github.com/temporalecologylab/bayes2025homework/blob/main/analyses/datasetsbayes.R). I checked with everyone and we have a good mix, so go forth with your planned dataset!

You can find some sample code for prior predictive checks [here](https://github.com/temporalecologylab/bayesclassdatasetsnmore/blob/main/analyses/misc/priorpc.R) and a complete (but simple) workflow for other data [here](https://github.com/lizzieinvancouver/bayesianflowsexample) (click on the example files). This latter workflow shows prior predictive checks and retrodictive checks. 

If you do retrodictive checks (aka posterior predictive checks), please challenge yourself to code one at least in R versus only looking in shinystan. 

At each step, please briefly describe what you did and why (e.g., for Step 1: Describe your model and how you decided on it), as well as noting any evaluations of your model that you did. Be sure to show you have examined your estimates and your model diagnostics for steps 2 and 4.  

### What file types can I submit?
I prefer .R or .Rmd for code, but Quarto works too. If you use R then you can write answers with comments (or submit a separate .md or .txt file about each step). Please include plots of your data and model (and perhaps both together!) and complete code. Compile your qmd and md files and submit the html or pdf with the .md and .qmd.

### What models should I fit?
Use rstanarm (`stan_glm` or `stan_lmer`). If your model does not converge **do not report the results** but do try to see if you can get it to run without divergences (check your code, consider an alternative model etc.). Rstanarm vignette [here](https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html).

## What if I need help?

Please connect with your peers in the class and/or post to [Piazza for the course](https://piazza.com/class/m5y7uqp7w133cj). I try to keep an eye on Piazza at least once a day (and, exceptionally, I will check it on the weekend once too). You're encouraged to help one another and work together, but please each submit a separate homework.  

# Due 27 January 2025 by 5pm

### Submit how? If you are enrolled in class, you must submit homework.
You can submit via Canvas or GitHub. To submit on GitHub:

1. Make a folder (name it your first name, all lowercase letters) in homeworksubmitted2/

2. Put your analyses and any output there, format the code so it runs from the repo without additional files.

3. Submit this as [a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request). 

If you submit on Canvas, follow steps 1-2, then compress your folder and upload it via the assignment. 

If you're auditing you don't have to submit your homework, but I still highly recommend doing it (it will test your skills in the course and make it easier for you to enjoy the homework review in class each week). You're welcome to submit your homework if you are auditing the course though please try to submit via one of the two methods above (I can add you to Canvas if you are auditing but need enough notice), because my email inbox is not a great place for homework. 

Thanks and have fun with this! 
