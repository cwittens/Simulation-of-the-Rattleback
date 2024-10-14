Hier sind alle Notebooks (Jupyter oder Mathematica) welche genutzt wurden um die Gleichungen bei Kuypers für die Winkelbeschleunigung (omega_dot), die Ableitung des Auflagepunktes (x_dot) zu berechnen.

Die Winkelbeschleunigung wurde in Jupyter wie auch in Mathematica berechnet und die daraus resultierende Gleichungen sind bis auf Rundungsfehler äquivalent. Die Gleichung für omega1_dot und omega2_dot in Mathematica war zwar kürzer (und damit auch schneller), aber leider nicht so stabil. Die Gleichungen für omega3_dot wurden von Mathematica genommen. Der Performance-Gewinn war aber nur minimal. 

Die beiden Mathematica Notebooks zum berechnen der naiven (insbesondere falschen!) Hamilton und Lagrange Gleichungen sind - to be best of my knowledge- korrekt. Da ich nicht denke, dass jemals irgendwer sich diese noch angucken wird, habe ich den Code nicht nochmal schön und übersichtlich gemacht - sorry. Falls doch und es fragen gibt, gerne bei mir melden :)


Die Lagrangefunktion aus *naive_lagrange_from_mathematica.nb* in Abhängigkeit von den Winkeln $\alpha,$ und $\beta$ und den Ableitungen $\dot{\alpha}, \dot{\beta}$ und $\dot{\gamma}$:
```math
ℒ = -\frac{g m \left(-h \cos (\alpha (t)) \cos (\beta (t)) \sqrt{\cos ^2(\alpha (t)) \left(a^2 \sin ^2(\beta (t))+c^2 \cos ^2(\beta (t))\right)+b^2 \sin ^2(\alpha (t))}+\cos ^2(\alpha (t)) \left(a^2 \sin ^2(\beta (t))+c^2 \cos
   ^2(\beta (t))\right)+b^2 \sin ^2(\alpha (t))\right)}{\sqrt{\cos ^2(\alpha (t)) \left(a^2 \sin ^2(\beta (t))+c^2 \cos ^2(\beta (t))\right)+b^2 \sin ^2(\alpha (t))}}+\frac{m \left(\left(b^2 \sin (\alpha (t)) \left(\dot{\alpha}(t)
   \sin (\beta (t))+\cos (\alpha (t)) \cos (\beta (t)) \dot{\gamma}(t)\right)-\left(\sin (\alpha (t)) \dot{\gamma}(t)+\dot{\beta}(t)\right) \left(c^2 \cos (\alpha (t)) \cos (\beta (t))-h \sqrt{\cos ^2(\alpha (t)) \left(a^2 \sin ^2(\beta
   (t))+c^2 \cos ^2(\beta (t))\right)+b^2 \sin ^2(\alpha (t))}\right)\right)^2+\left(\left(\dot{\alpha}(t) \cos (\beta (t))-\cos (\alpha (t)) \sin (\beta (t)) \dot{\gamma}(t)\right) \left(c^2 \cos (\alpha (t)) \cos (\beta (t))-h
   \sqrt{\cos ^2(\alpha (t)) \left(a^2 \sin ^2(\beta (t))+c^2 \cos ^2(\beta (t))\right)+b^2 \sin ^2(\alpha (t))}\right)+a^2 \cos (\alpha (t)) \sin (\beta (t)) \left(\dot{\alpha}(t) \sin (\beta (t))+\cos (\alpha (t)) \cos (\beta (t))
   \dot{\gamma}(t)\right)\right)^2+\left(\cos (\alpha (t)) \sin (\beta (t)) \left(\left(a^2-b^2\right) \sin (\alpha (t)) \dot{\gamma}(t)+a^2 \dot{\beta}(t)\right)+b^2 \dot{\alpha}(t) \sin (\alpha (t)) \cos (\beta (t))\right)^2\right)}{2
   \left(\cos ^2(\alpha (t)) \left(a^2 \sin ^2(\beta (t))+c^2 \cos ^2(\beta (t))\right)+b^2 \sin ^2(\alpha (t))\right)}+\frac{1}{2} \left(A \left(\dot{\alpha}(t) \cos (\beta (t))-\cos (\alpha (t)) \sin (\beta (t)) \gamma
   '(t)\right)^2+B \left(\sin (\alpha (t)) \dot{\gamma}(t)+\dot{\beta}(t)\right)^2+C \left(\dot{\alpha}(t) \sin (\beta (t))+\cos (\alpha (t)) \cos (\beta (t)) \dot{\gamma}(t)\right)^2+2 F \left(\sin (\alpha (t)) \dot{\gamma}(t)+\dot{\beta}(t)\right)
   \left(\cos (\alpha (t)) \sin (\beta (t)) \dot{\gamma}(t)-\dot{\alpha}(t) \cos (\beta (t))\right)\right)
```
