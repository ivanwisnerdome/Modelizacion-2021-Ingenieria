Hola, queria hacer un par de consultas del ejercicio 1 del TP 1.

Primero, trate de crear una carpeta para el TP 1, pero no supe como, en todo caso, subi el archivo.

Las preguntas son:

    Supongo que esta bien, compare los graficos que me dieron con los de la teoria, pero la verdad que no estaria muy seguro.

    En estos casos es para graficar los flujos, la verdad que no se supe muy bien como encararlo, y me copie de casi todo lo que hicieron en clase, osea lo que hizo Xavier. Trate de entenderlo todo, pero todavia no me quedo claro, por qué de las siguientes cosas:

¿ Por que la matriz de flujo[j,i] ? es decir porque los indices al reves

- porque el primer índice de la matriz corresponde a las filas, que se asocian 'hacia arriva' en el gráfico, a esa coordenada 
le pusimos j. la coordenada i es para la parte horizontal, como se acomodan las columnas. 

¿ Por que plt.streamplot(X,Y,-Flujox,-Flujoy)? ¿-flujox y -flujoy, de nuevo aca porque -flujox y -flujo en y

-de nuevo por lo mismo, depende también cómo hayas indexado tus coordenadas, pero la idea es la misma

    Es que creo que pude cambiar las condiciones de contorno, pero en la teoria claramente q=0 y el vector de carga tambien lo es, ¿ pero como cambia el vector de carga cuando el flujo es distinto de cero?, en lo que paso Mariano supongo que es bk=2deltaxqxA, no entiendo que es A aca.


Espero que se entienda, gracias, saludos

    tenes que ver como te cambia la ecuación cuando pones el flujo. fijate en la clase de ruben o en la diapo que subí al aula virtual, la ecuación para el punto ficticio:
$ \frac{\partialT}{\partial y} $= (algo proporcional a Q)
