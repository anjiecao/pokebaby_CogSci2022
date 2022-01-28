Welcome to our submission!

To reproduce our submission, please go to the 'writing' directory and open the 00_RANCH.Rmd file in RStudio. Then click 'Knit' at the top and it should render our submission as a pdf.

To run RANCH, along with the three different linking hypotheses (EIG, surprisal and KL),
please navigate to 'model' -> 'model_scripts' and open the 'model.Rproj' file to ensure you are running scripts from the right directory. From within RStudio, open minimal_simulation.Rmd and run each chunk. Note that to run the 500 simulations we report in our paper, especially for the EIG version of the model, will take a long time. So we recommend running it on a cluster, or running fewer simulations by adjusting n_sim in line 59 of the code.

To run our baseline models, please run the minimal_simulation_random_looking.Rmd and minimal_simulation_no_noise.Rmd scripts using the same R project window.
