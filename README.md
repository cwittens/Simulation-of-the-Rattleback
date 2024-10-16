# Simulation of the Rattleback
This is the GitHub Repository accompanying my mathematics Bachelor's thesis of simulating the rattleback. As the thesis was written in german, this README will also be in german from now on.

In diesem Repository findet sich der Code, um alle *Plots* in meiner Bachelorarbeit "Simulation des Wackelsteines" zu reproduzieren, die *Plots* selbst und noch weitere *Plots*, welche es nicht in die Bachelorarbeit selbst geschafft haben. Außerdem hat man die Möglichkeit verschiedenste Parameter, wie etwa die Anfangsbedingungen, die Toleranzen oder die Methoden zum Lösen der Bewegungsgleichen zu ändern. 

Der Quellcode der Verschiedenen Modelle findet sich in den Datein *kuypers.jl*, *schoemer.jl*, *naive_lagrangian.jl* und *naive_hamilton.jl*.

Um die *Plots* eines gewissen Kapitels zu reproduzieren, reicht es, das komplette Skript des Kapitels einmal unverändert auszuführen. Zum Beispeil in dem man im Julia REPL ``include("Kapitelnummer_Kapitelname.jl")`` ausführt. 

Um also Beispielsweise die *Plots* aus dem Kapitel "5.2 Stabilität der Modelle" zu reproduzieren, muss man im Julia REPL nur ``include("5.2_Stabilität der Modelle.jl")`` ausführen.

## Ein etwas anderer Abstract dieser Bachelorarbeit

O du Geist der Wahrheit, der durch die stürmischen Zeiten mich beseelt! Vernehmet nun die Kunde von einem wundersamen Steine, dem keltischen Wackelstein, dem Rattleback, jener eigentümlichen Erscheinung, die das stille Auge des Forschers fesselt und den Sinn des Denkers in schwindelnde Höhen erhebt.

Die vorliegende Arbeit, die dem ehrwürdigen Felde der Mathematik gewidmet ist, strebt danach, jenes seltsame Geschöpf der Natur, den Wackelstein, in der Sprache der Zahlen und Formen zu ergründen. Ein Halbellipsoid ist er, doch nicht eines von bloßer Vollkommenheit, nein, er birgt in seinem Innern eine geheimnisvolle Unordnung, eine inhomogene Verteilung seiner Massen, die ihn befähigt, die gewohnten Gesetze der Bewegung auf geheimnisvolle Weise zu durchbrechen und seine Drehrichtung plötzlich umzukehren. Ein Rätsel, das sich nicht auf bloßes Augenscheinliches zurückführen lässt, doch den Lösungen seiner erhabenen Bewegungsgleichungen gehorsam folgt.

<img src="https://upload.wikimedia.org/wikipedia/commons/f/f5/Celt_with_weights_of_gemstone_turtles-01.jpg" alt="Celt with weights of gemstone turtles" width="500"/>

Ein Bild des Wackelsteines, jenes widerspenstigen Geschöpfes! Die unheilvolle Verteilung der Massen wird hier durch die zusätzlichen Gewichte der Schildkröten bewirkt.

Es sei uns gestattet, mit bescheidenem Geiste einzutreten in die Hallen der Mechanik, wo wir die ehrwürdigen Gesetze, die uns Newton vermachte, auffrischen wollen. 

Wie uns Elmar Schömer lehrt, blicken wir zuerst auf diesen Stein durch das kalte Auge des Labors, wo er, gefangen im strengen Netz von 13 unerbittlich gekoppelten Differentialgleichungen, seine Bahn im Raum nach den eisernen Gesetzen der Mechanik offenbart.

Doch vernehmet auch die Stimme des weisen Friedhelm Kuypers, der uns den Wackelstein von einer anderen Warte aus zeigt, wie aus dem Inneren des Körpers selbst, und so enthüllt sich ein weiteres System, kleiner in der Zahl, doch nicht minder erhaben – sechs Gleichungen sind es, die die Gesetze seiner innersten Kräfte schildern.

So trete der Verstand ein in das Labyrinth dieser beiden Systeme! Wir wollen sie vergleichen, sie auf gleiche Pfade führen und prüfen, ob sie uns nicht gleichwohl zum selben Ziele tragen.

Denn gleichwohl sich die Form unterscheidet, bleibt die Wahrheit unveränderlich und fordert uns, die beiden Sichtweisen gegeneinander abzuwägen, auf dass wir den Stein nicht nur tiefer verstehen, sondern auch die Wege der Optimierung finden mögen, die uns zu neuen Höhen der Erkenntnis führen.

Und so, am Ende dieser Reise, wollen wir mit klarem Blick prüfen, welches der beiden Modelle die Gesetze der Energieerhaltung treuer wahrt, welches in den numerischen Lösungen beständiger bleibt und uns mit größerer Effizienz den Schlüssel zu den geheimnisvollen Tiefen des Wackelsteines in die Hand legt.