---
title: "Introduction"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---
FIXME

{% include links.md %}




Expectations- how big is my genome (and so how much sequence data should I generate?). 
    NB. For long reads- library properties have some influence on this.
    
    Also, cleaning up input data- how much of this do you want to do "live" vs have pre-prepared. Illumina I have taken through TrimGalore and then also got rid of everything >150bp since this makes the kmer profiling neater. 
    For ONT, I think I'd start with reads from which adaptors and control strands have already been removed, mostly to speed things up and focus time on the assembly. 