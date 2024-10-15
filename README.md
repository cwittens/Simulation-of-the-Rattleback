# Simulation of the Rattleback
This is the GitHub Repository accompanying my mathematics Bachelor's thesis of simulating the rattleback. As the thesis was written in german, this README will also be in german from now on.

In diesem Repository findet sich der Code, um alle *Plots* in meiner Bachelorarbeit "Simulation des Wackelsteines" zu reproduzieren, die *Plots* selbst und noch weitere *Plots*, welche es nicht in die Bachelorarbeit selbst geschafft haben. Außerdem hat man die Möglichkeit verschiedenste Parameter, wie etwa die Anfangsbedingungen, die Toleranzen oder die Methoden zum Lösen der Bewegungsgleichen zu ändern. 

Der Quellcode der Verschiedenen Modelle findet sich in den Datein *kuypers.jl*, *schoemer.jl*, *naive_lagrangian.jl* und *naive_hamilton.jl*.

Um die *Plots* eines gewissen Kapitels zu reproduzieren, reicht es, das komplette Skript des Kapitels einmal unverändert auszuführen. Zum Beispeil in dem man im Julia REPL ``include("Kapitelnummer_Kapitelname.jl")`` ausführt. 

Um also Beispielsweise die *Plots* aus dem Kapitel "5.2 Stabilität der Modelle" zu reproduzieren, muss man im Julia REPL nur ``include("5.2_Stabilität der Modelle.jl")`` ausführen.
