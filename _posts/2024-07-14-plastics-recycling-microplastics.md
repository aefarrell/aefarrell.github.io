---
title: "Plastics Recycling and Microplastics"
toc: true
toc_label: "Contents"
toc_sticky: true
comments: true
categories:
  - notes
tags:
  - plastic
  - recycling
  - microplastic
tagline: "is plastic recycling a huge source of microplastics?"
header:
  overlay_image: /images/plastics-recycling_files/soren-funk-unsplash.jpg
  overlay_filter: 0.25
  caption: "Photo by [Sören Funk](https://unsplash.com/@soerenfunk?utm_content=creditCopyText&utm_medium=referral&utm_source=unsplash) on [Unsplash](https://unsplash.com/photos/green-and-brown-stone-on-brown-sand-jQuky0VINAI?utm_content=creditCopyText&utm_medium=referral&utm_source=unsplash)>"
---

# Plastics Recycling and Microplastics: Checking the Math

As perhaps just a hazard of my profession, any time an article comes out on the merits (or lack of) of recycling and plastic waste in general, people send it my way. Several times in the last month I was sent [this article in *Quillette*](https://quillette.com/2024/06/17/recycling-plastic-is-a-dangerous-waste-of-time-microplastics-health/)<a href="#fn-1" class="sidenote-number"></a><span class="sidenote" id="fn-1">[Celia](#celia-2024), "Recycling Plastic is a Dangerous Waste of Time."</span>. (and associated [YouTube](https://www.youtube.com/watch?v=cJpkRSSM_OQ) video) about how plastics recycling may be a *massive* source of the microplastics being discharged into the environment, adding to the long list of reasons why recycling has not lived up to the promises made by industry, and undermining our path towards a more circular economy. At first glance though, some of the numbers presented and the math struck me as rather sus, so I would like to take a moment to dive into it a bit more. tl;dr much the math in that essay doesn't really work or comes with big caveats, but the broader point about the value of recycling and how we may not be fully appreciating the environmental impacts may hold up.

I won't be commenting on the broader life-cycle of plastic as I'm hardly an impartial participant: My current employer is one of the world's largest plastic manufacturers and my day job is working on a major expansion of one of its plastics manufacturing facilities. My employer has a whole suite of messaging about the importance of recycling and goals for advancing a circular economy which I don't feel particularly compelled to try and advance, nor contradict, on my little side project blog about doing math in my spare time. 
{: .notice }

## How Large of a Source of Microplastics is the Recycling Industry?

Celia estimates that up to 2/3<sup>rds</sup> of the microplastics discharged directly into the environment<a href="#fn-2" class="sidenote-number"></a><span class="sidenote" id="fn-2">[Celia](#celia-2024), "Recycling Plastic is a Dangerous Waste of Time." <br />These are so called "primary" microplastics, as opposed to "secondary" microplastics which are generated from plastic waste already in the environment</span> come from the recycling industry. This is a huge number. One that should immediately raise eye-brows. So lets break that down, it comes from two numbers:
1. that the recycling industry discharges up to 2Mt/y of microplastics into the environment
2. that the total amount of primary microplastics discharged into the environment from all sources is 3Mt/y

### The Direct Discharge of Microplastics

Celia takes a value of about 2Mt/y of microplastics from the recycling industry from an interview given by one of the authors of the study and leaves it rather mysterious as to where exactly that number comes from. However, this is a really easy number to calculate yourself: Approximately 9% of total plastic waste, globally, is recycled, elsewhere in the essay Celia claims 6 - 13% of recycled plastic is lost to the environment as primary microplastics, which equates to about 2 - 4Mt/y.

For example, the OECD estimated that, in 2019, global plastic waste generation was 353Mt of which 33Mt were recycled (~9%)<a href="#fn-3" class="sidenote-number"></a><span class="sidenote" id="fn-3">[OECD](#oecd-2022), *Global Plastics Outlook*, 20.</span> 6 - 13% of that is 1.98 - 4.29 Mt. So in some sense, taking the high end of that, makes the argument *more dramatic*.

The main reference around which the entire essay revolves is a recent study of a single plastics recycling plant in the UK<a href="#fn-4" class="sidenote-number"></a><span class="sidenote" id="fn-4">[Brown *et al*](#brown-2023), "Potential for a Plastic Recycling Facility to Release Microplastic Pollution."</span>. In this study, the authors looked at a relatively new plastics recycling facility that underwent an upgrade to its wastewater treatment process, adding additional filtration. The study looked at the microplastics prior to and after that upgrade. Based on the concentrations measured in the wastewater they estimated that up to 13% of the mass of plastic brought into the facility may have been discharged in the wastewater as microplastics prior to the filtration upgrade and, after the upgrade, this dropped to 6%. These two numbers 6% and 13% form the basis for the estimate of how much primary microplastics are being discharged from the recycling industry as a whole.

At this point we should pause consider *where those numbers came from*. The study gives a range for the total annual mass discharge in the wastewater, based on measured mass concentrations in the facility wastewater. The ratio of this mass out to the plastic taken in is the origin of the 6 - 13% range. However, I think it is deeply disingenuous to present these numbers without context as the study's estimates span *three orders of magnitude*.

| estimate              | low end (t/y) | high end (t/y) |
|:----------------------|:-------------:|:--------------:|
| before filter upgrade | 96            | 2933           |
| after filter upgrade  | 4             | 1366           |

I think the take away from this is that far more data is needed to narrow these error bars. The low end estimates are in line with other studies for the whole life-cycle of plastic<a href="#fn-5" class="sidenote-number"></a><span class="sidenote" id="fn-5">[Ryberg *et al*](#ryberg-2018), *Mapping of Global Plastics Value Chain*, 56; for example assumes that microplastic losses are ~0.005% for recycling, which is in the same range as the low end estimate for the facility prior to the filtration upgrade</span>, the high end estimates are many orders of magnitude larger. This isn't to say that these values are wrong, merely that it reveals *the potential* for far more microplastic losses to the environment from the recycling industry than have previously been accounted for.


### The Total Amount of Primary Microplastics losses

I think the ~3Mt/y is a relatively robust estimate, for the type of study Celia references, because it has been replicated<a href="#fn-6" class="sidenote-number"></a><span class="sidenote" id="fn-6">[Ryberg *et al*](#ryberg-2018), *Mapping of Global Plastics Value Chain*, estimates 3.1Mt/y;<br /> [OECD](#oecd-2022), *Global Plastics Outlook*, estimates 2.7Mt/y;<br /> [Boucher and Friot](#boucher-2017), *Primary Microplastics in the Oceans*, estimates 1.8 - 5Mt/y.</span>. However, this is the source of the most egregious and obvious mistake, and the one that prompted me to write this blog post in the first place. The studies referenced as the sources for the 3Mt/y number include recycling as a source in the estimate but *do not estimate the losses from recycling to be anywhere near as high*. Dividing these two numbers is entirely mathematically invalid.

Before I go any further, where do numbers like 3Mt/y come from? They are not from direct measurements of microplastics in the environment. They come from a life-cycle analysis that looks at the entire life of plastics and estimates rates of losses at the various steps along the path from the creation of virgin plastic to its ultimate fate. Adding all of these losses up gives the total estimated primary microplastics loss. This is why it is incorrect to simply ratio 2Mt/y over 3Mt/y: that would only work if the 2Mt/y was *included in the total*, and it isn't.

Supposing that we are going with 2Mt/y of primary microplastics from recycling, most studies (importantly the ones referenced by Celia) do not use a number anywhere near this high. In fact most assume it quite small and some take it to be negligible<a href="#fn-7" class="sidenote-number"></a><span class="sidenote" id="fn-7">[Ryberg *et al*](#ryberg-2018), *Mapping of Global Plastics Value Chain*, 56.</span>. The correct procedure would be to *subtract* the previous estimate for losses due to recycling from the total losses, and then add the new estimate of 2Mt/y, giving a corrected total. This would then be the denominator.

Consider a UNEP study that estimates the total primary microplastics losses from the entire plastics value chain as 3.1Mt/y<a href="#fn-8" class="sidenote-number"></a><span class="sidenote" id="fn-8">[Ryberg *et al*](#ryberg-2018), *Mapping of Global Plastics Value Chain*, 54.</span> Conveniently, this study assumed the losses due to recycling were negligible (i.e. zero). So based on this study's estimates for *all other sources of primary microplastics*, and our estimate of 2Mt/y from the recycling industry alone, we would estimate a *new total* of 5.1Mt/y, of which 2/5.1 = 39% came from the recycling industry. So.. not 2/3<sup>rds</sup>.

But adding in how wide the error bars are for the estimate of how much primary microplastics are discharged into the environment, from that one plastics recycling facility, all we can really say is that recycling is somewhere between a negligible source of primary microplastics and the single largest source of primary microplastics. Which is important in the sense that it identifies that we *may* be missing a major source of primary microplastics, but it really does not live up to the hype in Celia's article.

## The Climate Impacts of Landfilling

Celia makes reference to the landfilling of municipal solid waste being a source of methane emissions as part of the argument for why recycling should be abandoned and plastic incinerated instead. Independent of the merits of recycling or incinerating, this is at best irrelevant. Plastic has a negligible methane generating potential when landfilled, a fact that is related to the primary concern with plastic waste in the environment: its environmental persistence. The methane emissions coming from the landfilling of municipal waste is from decomposing organic matter, not the plastic. In fact a recent meta-analysis<a href="#fn-9" class="sidenote-number"></a><span class="sidenote" id="fn-9">[Gao *et al*](#gao-2024), "Comprehensive Meta-Analysis."</span> shows that, if anything, the presence of non-biodegradable plastic *reduces* the methane emissions from anaerobic digestion as non-biodegradable plastics may leach toxins that prevent bacteria from decomposing organic matter. I wouldn't take that to mean we *should be* landfilling plastic waste, as some climate mitigation strategy, merely that the methane emissions from doing so are irrelevant to the argument around what to do with plastic waste.

## Conclusions and Take Aways

There is certainly a growing chorus of concern over the fate of plastics in the environment, and the environmental and health consequences of microplastics given their ubiquity. That alone should warrant a lot more study into the sources of microplastics. That the estimate that recycling accounts for 2/3rds of primary microplastics doesn't hold up, due to rudimentary math mistakes, doesn't invalidate the broader concern that recycling simply has not lived up to the promise and may in fact be worsening the microplastics problem. We don't know that is the case, given the data cited, but I think the onus is on the recycling industry to show that they are, in fact, part of the solution and not making the problem worse.

I am not going to comment on the relative merits of incineration, recycling, or advanced recycling other than to say few of the technical problems in this field are truly insurmountable. The real question always comes down to *cost* and how much we are willing to pay to achieve the environmental performance we want.


## References

+ <a name="boucher-2017">Boucher</a>, Julien, Damien Friot. *Primary Microplastics in the Oceans*. Gland, Switzerland: International Union for Conservation of Nature. 2017. [doi: 10.2305/IUCN.CH.2017.01.en](https://doi.org/10.2305/IUCN.CH.2017.01.en)
+ <a name="brown-2023">Brown</a>, Erina, Anna MacDonald, Steve Allen, Deonie Allen. "The Potential for a Plastic Recycling Facility to Release Microplastic Pollution and Possible Filtration Remediation Effectiveness." *Journal of Hazardous Materials Advances*. 10 (2023): 100309. [doi: 10.1016/j.hazadv.2023.100309](https://doi.org/10.1016/j.hazadv.2023.100309)
+ <a name="celia-2024">Celia</a>, Frank. "Recycling Plastic is a Dangerous Waste of Time." *Quillette*. June 17, 2024. [url](https://quillette.com/2024/06/17/recycling-plastic-is-a-dangerous-waste-of-time-microplastics-health/)
+ <a name="gao-2024">Gao</a>, Zhenghui, Hang Qian, Tianyi Cui, Zongqiang Ren, Xingjie Wang. "Comprehensive Meta-Analysis Reveals the Impact of Non-Biodegradable Plastic Pollution on Methane Production in Anaerobic Digestion." *Chemical Engineering Journal*. 484 (2024): 149703. [doi: 10.1016/j.cej.2024.149703](https://doi.org/10.1016/j.cej.2024.149703)
+ <a name="oecd-2022">OECD</a>. *Global Plastics Outlook*. Paris: OECD Publishing. 2022. [doi: 10.1787/de747aef-en](https://doi.org/10.1787/de747aef-en)
+ <a name="ryberg-2018">Ryberg</a>, Morten W., Alexis Laurent, Michael Hauschild. *Mapping of Global Plastics Value Chain and Plastics Losses to the Environment*. Nairobi: United Nations Environment Programme. 2018. [url](https://www.unep.org/resources/report/mapping-global-plastics-value-chain-and-plastics-losses-environment-particular/)
