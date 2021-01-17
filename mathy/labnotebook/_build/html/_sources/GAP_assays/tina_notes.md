# Notes from Tina about GAP assays

- Tina would do one replicate with/without GAP (so two wells total) per loading, "Because that was the largest variation in values. Experiments from the same fractions of loading set up in the same way were really the same."
    - She says 2 technical replicates for GAP and 2 for no GAP per loading would be good, but sensor is the biggest expense so she just did one.
    - 3 "biological" replicates which were the GTP loadings

- If you don't know loading efficiency, you may want to do 8 uM or 10 uM so that final concentration is closer to 5 uM
    - "it seems lower conc were better - lower fitting error"
    - "So 5 is good but also use 10 if your loading efficiency is low so you get 5 in the end"

- prewarming the plate is important for the start being linear

- use 20 uM sensor
    - "I used to play with the ratio of sensor and protein. I guess with 10 we could have used 20 not 50 of sensor"
    - "You want to be up to 60% of sensor" i.e. 0.6*[sensor conc.] = [Gsp1 conc]
    - "If your phosphate concentration is over 60% sensor concentration you are out of the linear range"

- on the calibration:
    - "So the calibration I have is linear. If we order new sensor we should order some more and calibrate the 20 uM. Ideally you want to do everything with the same lot and have that lot calibrated. I think this is where most of our error came from. Mixing lots and calibrations.
    - "Also you don’t need to calibrate before you run. You can calibrate any time. It’s to get final concentration. And you need that to get km and kcat. But for first approx you can use my old calibration values."
    - I asked "i do need to re-calibrate for every lot, right?", to which she responded "Yes ideally. I think I had two of three calibrated in the end. Cos it uses up so much sensor. And the values were pretty similar. But I think it’s a good idea to calibrate and hopefully you will do everything with the same lot. I think the main benefit now is that we learned it’s ok to do the experiment with low concentration of loaded protein. Like 8. So you can do all the exp with 20 uM sensor. I had a bunch of calibrations for 10, 20, 50, 70 uM. And that ate up multiple tubes of the sensor. You don’t have to do that"

- n. of experiments with 3 tubes of sensor:
    - 200 ul per tube * 500 uM
    - dilute it 4 times to 125 uM that makes 800 ul
    - you use 20 of that per reaction
    - So that is 40 reaction per tube
    - You want to do 30 reactions at least plus calibration
    - we can easily end up using 3 tubes

- Tina after fitting the F28V and F28Y curves
    - The curves were a bit bumpy, especially the 10 uM F28V was hard to fit, so I would't trust that super low Km. The rest seems OK. looks like the tolerant mutant is very similar to WT and the DN mutant is just a bit slower. in general not a large effect!

- On minimizing batch variation:
    - I said: "Hey can I run something by you regarding GAP assay? So I don’t have as many aliquots frozen as we did when we had larger scale preps, so Tanja and I figured doing two loadings in parallel from one thawed tube can count as a “biological replicate” for our purposes. Do you agree? So I would thaw one tube, do two loadings side by side, then with each loading do +/- GAP."
    - She said: another thing you can add to make it more viariable, or more like it's two different days, is use two separate sensor dilutions for the two batches. I think you will get more variations from that than from loading, to be fair. so lets say you do samples 1-3: you run 1a, 2a, 3a with one dilution of sensor and 1b, 2b, and 3b with another 4x dilution of sensor. do you see what I mean?"
    - I said: "yes, i do, so the suggestion is that most of our variability across the entire dataset will be from different sensor dilutions". She said "that's my feeling". and I said "so to avoid a day-based batch effect, I should use two sensor dilutions on the same day, but i should still do two loadings, since that will account for variability between different loading efficiencies.

- On loading times:
    - Tina also said I could leave one of the two loadings on the NAP5 column for longer, to simulate different loading times. I mentioned I'm less worried because it's essentially at equilibrium after 3 hours, given minimal hydrolysis on ice. She agreed but mention it could also take aggregation into account.
    - I agreed that I could have 6 NAP5s for 3 samples x 2 loadings, and elute 1a, 2a, 3a first, do my nanodrop measurements which takes some time, then elute 1b, 2b, 3b.

- Sensor dilutions:
    - I asked: "should i aliquot it today when i thaw it? like aliquot into black tubes, then i can just take a tube out and immediately dilute with a certain volume?"
    - Tina's response: "yes, I used to aliquot after thawing the first time because freeze thawing definitely affects it. I don't dilute it before freezing. always flash freeze in it's original buffer, and at high conc. I always dilute it with the same buffer I use for assay. I would write down on the side of the black tube the volume, like 32 ul and then once I took it out, I would just spin down and add buffer on top of it. so make sure to calculate how many runs you want to use and how big your aliquot needs to be.
    - She then proceeded: "the most important thing is to use the same diluted sensor for +/- GAP. diluting the sensor twice 4x is a huge variation, you will see. don't mix and match sensor for +/- GAP". I clarified that you can tell from the starting point of the curves being the same, and she confirmed.