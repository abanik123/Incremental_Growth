# Paper Title: Modelling incremental length-based growth using a Cormack–Jolly–Seber framework: An application to anadromous Dolly Varden
<h1>Project Overview</h1>

<p>
This repository contains R code, simulation scripts, and datasets used for modeling the growth dynamics 
of Dolly Varden under multiple biological scenarios. The project is organized into three main components:
</p>

<ul>
    <li><strong>Case_Study</strong> – R code for fitted models applied to two river systems.</li>
    <li><strong>Simulation_Study</strong> – Simulation scripts for evaluating model performance.</li>
    <li><strong>Data</strong> – Datasets corresponding to the case studies.</li>
</ul>

<p>
The models implemented include:
</p>
<ul>
    <li>Unconditioned growth</li>
    <li>Growth depending on breeding state</li>
    <li>Growth depending on sex</li>
</ul>

<hr>

<h2>1. Case_Study Codes/</h2>

<p>This directory contains analyses for two river systems:</p>

<ul>
    <li><strong>BigFish/</strong></li>
    <li><strong>JoeCreek/</strong></li>
</ul>

<p>Each river folder includes the following subfolders:</p>

<ul>
    <li><strong>SAGA/</strong></li>
    <li><strong>VBGM/</strong></li>
</ul>

<p>
Both SAGA and VBGM contain R code implementing the three growth modeling cases:
</p>

<ul>
    <li>Unconditioned growth</li>
    <li>Growth depending on breeding</li>
    <li>Growth depending on sex</li>
</ul>

<p>
These scripts reproduce the case-study results presented in the associated research paper.
</p>

<hr>

<h2>2. Simulation_Study Codes/</h2>

<p>
This directory includes R scripts used to conduct simulation studies.
The provided scripts simulate a baseline scenario; however, users can modify the parameter values 
to match other scenarios using the values described in the paper.
</p>

<hr>

<h2>3. Case_Study Data/</h2>

<p>
This folder contains datasets used in the case studies and includes:
</p>

<ul>
    <li><strong>Big Fish River/</strong></li>
    <li><strong>Joe Creek/</strong></li>
</ul>

<p>
Each river folder contains datasets corresponding to the three model cases:
</p>

<ul>
    <li>Unconditioned growth</li>
    <li>Growth depending on breeding</li>
    <li>Growth depending on sex</li>
</ul>

<hr>

<h2>Requirements</h2>

<ul>
    <li>R (4.4.0)</li>
    <li>R packages commonly used for Bayesian and growth modeling (e.g., <code>nimble</code>, <code>MCMCvis</code>, <code>tidyverse</code>,
      <code>dplyr</code>, <code>abind</code>)</li>
</ul>

<hr>

<h2>Usage</h2>

<ol>
    <li>Navigate to <strong>Case_Study/</strong> to reproduce analysis for Big Fish River or Joe Creek.</li>
    <li>Navigate to <strong>Simulation_Study/</strong> to run simulation experiments or modify parameters.</li>
    <li>Use the datasets in <strong>Data/</strong> when running or adapting scripts.</li>
</ol>

<hr>
