This file last update on 19/03/19

This protocol is designed to detect events for each species separately.
Serves a secondary role as reference for the best species-specific criteria (including the two levels of click
discrimination plus event detection).

Developed this in light of the realization that event detection criteria may have very different sensitivities 
for different species (see investigation notes below).

Criteria sheets with highlighting note which criteria and thresholds are different from the general set.


CAVEATS:
I've only looked at two sets of event detection criteria. For those where 5% beaked was shown to work, there
could always be better systems that improve precision.
- "Standard" set (>= 5% beaked)
- "Basic" set (>= 1 beaked)


FINDINGS ON EVENT DETECTION:
- Event Detection criteria for Mb work much better with simplest case (i.e. need only 1 click). However, for Ha, 
5% beaked is better.
- Perhaps this can be justified: Mb clicks are more unique (less overlap with dolphins), so need fewer of them to
be convinced.
- 5% beaked was not an issue for Mb with general criteria. Investigation shows that this occurs because of Ha.
Mb is often found with Ha (82% in my Midgul data so far), so the detector with 5% beaked was mostly picking up Ha,
and less so Mb. Logically, this poses a bit of a problem with the current "protocol" system, because in the case where I use multiple discriminators for an event, it kind of assumes that the same event detection criteria apply
to each set of clicks.
- Relaxing event-level criteria for Mb does not solve the problem. More evidence that Ha biased the Mb recall
greatly.
- It kind of makes sense actually. If detector doesn't look for Ha but the majority of clicks are Ha, Mb won't meet the 5% criterion.
- 5% beaked is fine for Me.


PROGRESS UPDATES:
Ha DECLARED FINAL ON 18/03/19
Mb DECLARED FINAL ON 18/03/19
Me DECLARED FINAL ON 19/03/19
Zc DECLARED FINAL ON 19/03/19