/*
Se pide desarrollar una solución, para armar palabras de ADN, así como su frecuencia de aparición; de
acuerdo con un tamaño k que deberá pedirle al usuario (valores de k deben de ser entre 4 y 10); una vez
armada las palabras se pide que se realice el conteo de aparición de cada palabra, así mismo se deberá
desarrollar un mecanismo de búsqueda de palabras, en el cual, el usuario proporcionará una palabra y el
programa deberá informar si la palabra se encuentra, cuál es su frecuencia de aparición.
Ejemplo
Usuario pide PALABRAS DE 5
Usando una cadena de ADN: AGCTTTTNCATTCTGACTGCAACGGGCAATATG (Se usarán archivo
fna o archivos fasta)
Se inicia de izquierda a derecha por la letra A y se toman 5 letras (tomando la primera de inicio) la
primera palabra seria AGTTT, la siguiente palabra se toma al correr la posición inicial en 1, es decir la
segunda palabra inicia en la posición 2 que en este caso es la G y se toman 5 letra quedaría GCTTTT, se
corre una posición y se toman 5 en este caso inicia con la C y quedaría la palabra CTTTT y así
sucesivamente hasta agotar la cadena del archivo.
La salida de las palabras sería algo así y sus frecuencias de aparición

Integrantes: 
Canche Ciau Rusell Emmanuel
Gutierrez Perez Claudio Habraham
*/