Real-Time Rendering of Glints in the Presence of Area Lights
============================================================

Overview
--------

This repository contains the source code to reproduce the figures from the article
> Tom Kneiphof and Reinhard Klein. "Real-Time Rendering of Glints in the Presence of Area Lights." Pacific Graphics 2024.

It contains a Unity3D project and a modified version of Unity3D's High Definition Render Pipeline containing our implementation [including parts of previous work](https://thomasdeliot.wixsite.com/blog/single-post/hpg23-real-time-rendering-of-glinty-appearance-using-distributed-binomial-laws-on-anisotropic-grids).
The Beethoven bust from the teaser figure can be found [here](https://www.cgtrader.com/3d-models/architectural/decoration/beethoven-bust-47f58cebea4e58aaa8946756726eb47d) and should replace `Assets/Models/Beethoven.obj`.


Project Structure
-----------------

A copy of the HDRP and core render pipelines are found in the `Packages/` directory.
Our implementation modifies the `HDRP/Lit` shader in the HDRP.
The relevant entry points for the glint evaluation are found in the functions `DirectLighting EvaluateBSDF_Rect(...)` and `CBSDF EvaluateBSDF(...)` in `Packages/com.unity.render-pipelines.high-definition@14.0.10/Runtime/Material/Lit/Lit.hlsl`.
The actual implementation of our method and previous work is split across the files `Packages/com.unity.render-pipelines.high-definition@14.0.10/Runtime/Material/Glints/Glints{2023,2024,Subdivision}.hlsl`.

The Unity Project itself contains the scenes to reproduce figures and some utility scripts.


Reproducing Figures
-------------------

By default, the scenes in `Assets/Scenes/.../` produce the raw screenshots that we used to assemble our figures.
Each scene contains a "TriggerScreenshots" object, responsible for initiating the capturing of screenshots, that will be placed in the `Output/` directory.

An exception is the "RepresentativeVideo" scene, for which video caputre is achieved using the UnityRecorder in a Timeline after starting the play mode.
Individiaul frames will be placed in `Recordings/` where a Makefile can assemble the frames into a video.

Another exception is the "RuntimePlane" scene for which an executable should be built that creates csv files with runtime measurements upon execution.


| Scene                | Figures in Paper |
|----------------------|------------------|
| AreaLightRefAblation | Figs. 7, 8       |
| AreaLightSubdivision | Fig. 9           |
| Beethoven            | Fig. 1           |
| Convergence          | Fig. 10          |
| MicrofacetRoughness  | Fig. 6           |
| RepresentativeVideo  | Supp. Video      |
| RuntimePlane         | Tab. 1           |


Citation
--------

If you are using this code in academic research, please cite our paper.
The BibTeX entry is (subject to change)
```bibtex
@article{Kneiphof2024Real,
    author = {Tom Kneiphof and Reinhard Klein},
    title = {Real-Time Rendering of Glints in the Presence of Area Lights},
    journal = {Computer Graphics Forum},
    volume = {43},
    number = {7},
    month = {oct},
    doi = {TBD}
}
```
