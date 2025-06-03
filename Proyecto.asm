; Proyecto.asm
; Programa para analizar palabras de ADN en archivos FASTA/FNA
; Características:
;   - Lee archivos FASTA/FNA, omite líneas de cabecera
;   - Genera palabras de longitud k (4-10)
;   - Guarda SOLO palabras ordenadas con repeticiones en "palabras.txt"

section .data                         ; Sección de datos inicializados en memoria

LF equ 10                             ; Constante: salto de línea (Line Feed) ASCII = 10
NULL equ 0                            ; Constante: valor nulo (terminador de cadena)
EXIT_SUCCESS equ 0                    ; Código de salida para éxito (0)
STDIN equ 0                           ; Descriptor de archivo para entrada estándar (teclado)
STDOUT equ 1                          ; Descriptor de archivo para salida estándar (pantalla)
STDERR equ 2                          ; Descriptor de archivo para salida de errores
SYS_read equ 0                        ; Número de syscall para leer desde archivo/dispositivo
SYS_write equ 1                       ; Número de syscall para escribir a archivo/dispositivo
SYS_open equ 2                        ; Número de syscall para abrir un archivo
SYS_close equ 3                       ; Número de syscall para cerrar archivo
SYS_exit equ 60                       ; Número de syscall para terminar el proceso
SYS_creat equ 85                      ; Número de syscall para crear archivo
O_RDONLY equ 000000q                  ; Abrir archivo solo lectura
O_WRONLY equ 000001q                  ; Abrir archivo solo escritura
O_CREAT equ 0x40                      ; Crear archivo si no existe
O_TRUNC equ 0x200                     ; Truncar archivo si ya existe
O_APPEND equ 0x400                    ; Escribir al final del archivo
S_IRUSR equ 00400q                    ; Permiso: lectura para el usuario
S_IWUSR equ 00200q                    ; Permiso: escritura para el usuario
BUFF_SIZE equ 65536                   ; Tamaño del buffer de lectura = 64KB (2^16)
MAX_ADN_SIZE equ 8000000              ; Tamaño máximo de ADN que se almacenará = 8,000,000 bytes (8MB)  
separator db ": ", 0                  ; Cadena que separa k-mer y frecuencia en la salida (": ")
header db "Analizador de palabras de ADN", LF, LF, NULL ; Mensaje de encabezado que se muestra al iniciar
fileName db "Prueba.txt", NULL        ; Nombre del archivo de entrada
outFileName db "palabras.txt", NULL   ; Nombre del archivo de salida
promptK db "Ingrese longitud de palabra (4-10): ", NULL ; Mensaje que pide al usuario el valor de k
espacio db " ", 0                     ; Cadena con espacio simple 
ordenadoMsg db "Palabras Ordenadas con Repeticiones:", 13, 10, 0 ; Encabezado para la sección de k-mers ordenados
errMsgOpen db "Error abriendo archivo.", LF, NULL ; Mensaje de error si no se puede abrir el archivo de entrada
errMsgRead db "Error leyendo archivo.", LF, NULL ; Mensaje de error si falla la lectura del archivo
errMsgWrite db "Error escribiendo archivo.", LF, NULL ; Mensaje de error si falla la escritura al archivo
errMsgK db "Error: k debe ser entre 4 y 10.", LF, NULL ; Mensaje de error si el valor de k ingresado es inválido
successMsg db "Palabras guardadas en palabras.txt", LF, NULL ; Mensaje mostrado al finalizar exitosamente
crlf db 13, 10, 0                     ; Retorno de carro + salto de línea (CR + LF), usado al final de cada línea en archivo

section .bss                         ; Sección para datos no inicializados (se llenan con ceros en tiempo de ejecución)

readBuffer resb BUFF_SIZE            ; Buffer temporal para leer bloques del archivo (64KB)
adnBuffer resb MAX_ADN_SIZE          ; Buffer grande para almacenar toda la secuencia de ADN limpia (hasta 8MB)
fileDesc resq 1                      ; Variable para guardar el descriptor del archivo de entrada (8 bytes)
outFileDesc resq 1                   ; Variable para guardar el descriptor del archivo de salida (8 bytes)
kVal resb 1                          ; Valor de k (longitud de palabra a generar), ingresado por el usuario
adnLength resq 1                     ; Longitud total de la secuencia ADN cargada en adnBuffer
kmerList resb MAX_ADN_SIZE           ; Espacio para almacenar todos los k-mers generados, uno detrás de otro
kmerCounts resq MAX_ADN_SIZE         ; Contador de repeticiones por cada k-mer (paralelo a kmerList)
tempKmer resb 16                     ; Buffer temporal para comparar un k-mer durante el conteo (hasta 16 caracteres)
contador resq 1                      ; Variable auxiliar de uso general para ciclos o contadores
totalKmers resq 1                    ; Cantidad total de k-mers generados y guardados en kmerList


section .text
global _start

_start:
    mov rdi, header              ; Cargar en RDI la dirección del mensaje de cabecera "Analizador de palabras de ADN"
    call printStr               ; Llamar a la función que imprime la cadena por consola

    call openInputFile          ; Llamar a la función para abrir el archivo de entrada
    cmp rax, 0                  ; Verificar si hubo error (valor negativo)
    jl _exitError               ; Si error, salir del programa

    call readFASTA              ; Leer y procesar el archivo FASTA (filtrando solo caracteres válidos de ADN)
    cmp rax, 0                  ; Verificar si hubo error
    jl _exitError               ; Si error, salir del programa

    mov rax, SYS_close          ; Preparar syscall para cerrar archivo
    mov rdi, [fileDesc]         ; Obtener el descriptor del archivo abierto
    syscall                     ; Llamar al sistema para cerrar el archivo

    call getKValue              ; Solicitar al usuario la longitud de palabra (k)
    cmp rax, 0                  ; Verificar si hubo error de entrada
    jl _exitError               ; Si error, salir del programa

    call generateWords          ; Función "placeholder", puede validar o preparar datos
    cmp rax, 0                  ; Verificar si hubo error
    jl _exitError               ; Si error, salir del programa

    ; Crear archivo de salida (si ya existe, lo trunca a cero)
    mov rax, SYS_creat          ; Preparar syscall para crear archivo
    mov rdi, outFileName        ; Nombre del archivo de salida ("palabras.txt")
    mov rsi, S_IRUSR | S_IWUSR  ; Permisos: lectura y escritura para el usuario
    syscall                     ; Llamar al sistema para crear archivo
    cmp rax, 0                  ; Verificar si hubo error
    jl _exitError               ; Si error, salir del programa
    mov [outFileDesc], rax      ; Guardar el descriptor del archivo de salida

    call ordenarYGuardarKmers   ; Generar todos los k-mers a partir del ADN y guardarlos
    call ordenarKmers           ; Ordenar alfabéticamente los k-mers generados

    mov rdi, ordenadoMsg        ; Cargar mensaje "Palabras Ordenadas con Repeticiones:"
    call printStrToFile         ; Imprimirlo en el archivo de salida

    call contar_frecuencias     ; Contar y escribir la frecuencia de cada k-mer ordenado

    mov rax, SYS_close          ; Preparar syscall para cerrar archivo
    mov rdi, [outFileDesc]      ; Descriptor del archivo de salida
    syscall                     ; Cerrar archivo

    mov rdi, successMsg         ; Mensaje de éxito ("Palabras guardadas en palabras.txt")
    call printStr               ; Imprimirlo por consola

    mov rax, SYS_exit           ; Preparar syscall para terminar el programa
    mov rdi, EXIT_SUCCESS       ; Código de salida 0 (éxito)
    syscall                     ; Salir del programa

_exitError:                     ; Etiqueta en caso de error
    mov rax, SYS_exit           ; Preparar syscall de salida
    mov rdi, 1                  ; Código de salida 1 (error)
    syscall                     ; Salir del programa con error
;-------------------------------------------
; Función para imprimir cadenas
printStr:
    push rbx                 ; Guardamos el registro rbx en la pila (por si lo está usando quien llama)

    mov rbx, rdi             ; Guardamos el puntero a la cadena (que está en rdi) en rbx
    xor rdx, rdx             ; Inicializamos el contador de longitud en rdx a 0

.cont:
    cmp byte [rbx + rdx], 0  ; Comparamos el byte actual con 0 (fin de cadena)
    je .done                 ; Si es 0, terminamos el conteo
    inc rdx                  ; Si no, aumentamos el contador de longitud
    jmp .cont                ; Repetimos el proceso hasta encontrar el fin de cadena

.done:
    mov rax, SYS_write       ; Código de syscall para write (escribir)
    mov rdi, STDOUT          ; Descriptor del archivo de salida: STDOUT (pantalla)
    mov rsi, rbx             ; Dirección de inicio de la cadena a imprimir
    syscall                  ; Realizamos la llamada al sistema para imprimir

    pop rbx                  ; Restauramos el valor original de rbx
    ret                      ; Regresamos al llamador
;-------------------------------------------
; Función para imprimir cadenas en archivo
printStrToFile:
    push rbx                 ; Guardamos el registro rbx en la pila (backup)

    mov rbx, rdi             ; Copiamos el puntero a la cadena (desde rdi) en rbx
    xor rdx, rdx             ; Inicializamos el contador de longitud en 0

.lenloop:
    cmp byte [rbx + rdx], 0  ; Comparamos el byte actual con 0 (NULL terminador)
    je .write                ; Si encontramos el fin de la cadena, saltamos a escribir
    inc rdx                  ; Si no, aumentamos el contador de longitud
    jmp .lenloop             ; Repetimos hasta encontrar el NULL

.write:
    mov rax, SYS_write       ; Código de syscall para escribir (write)
    mov rdi, [outFileDesc]   ; Cargamos el descriptor del archivo de salida (palabras.txt)
    mov rsi, rbx             ; Dirección de la cadena a escribir
    ; rdx ya tiene la longitud desde el bucle anterior
    syscall                  ; Llamamos al sistema para escribir la cadena en el archivo

    pop rbx                  ; Restauramos el valor original de rbx
    ret                      ; Regresamos al llamador
;-------------------------------------------------------
; Copiar k-mers en kmerList para ordenarlos
ordenarYGuardarKmers:
    xor rsi, rsi                 ; rsi = 0, índice para recorrer adnBuffer (inicio de lectura)
    xor rdi, rdi                 ; rdi = 0, índice para guardar en kmerList

    movzx rcx, byte [kVal]       ; rcx = longitud del k-mer (valor ingresado por el usuario)
    mov rbx, [adnLength]         ; rbx = longitud total del ADN leído

    cmp rbx, rcx                 ; ¿Hay suficientes caracteres para al menos un k-mer?
    jb .fin                      ; Si no, salir (no hay suficientes datos)

    sub rbx, rcx                 ; rbx = total posibles posiciones para k-mers
    inc rbx                      ; Agregamos 1 porque se puede formar un k-mer en la última posición válida
    mov [totalKmers], rbx        ; Guardamos cuántos k-mers vamos a generar

.copiar:
    cmp rbx, 0                   ; ¿Ya copiamos todos los k-mers?
    je .fin                      ; Si sí, terminamos

    mov rdx, 0                   ; rdx será el índice para copiar cada letra del k-mer

.loop_kmer:
    cmp rdx, rcx                 ; ¿Ya copiamos k letras?
    je .next                     ; Si sí, pasamos al siguiente k-mer

    mov al, [adnBuffer + rsi + rdx]     ; al = letra de ADN actual
    mov [kmerList + rdi + rdx], al      ; Copiamos la letra al arreglo de k-mers
    inc rdx                              ; Avanzamos al siguiente carácter del k-mer
    jmp .loop_kmer                       ; Repetimos hasta tener el k-mer completo

.next:
    inc rsi                      ; Avanzamos una posición en adnBuffer (para siguiente k-mer)
    add rdi, rcx                 ; Avanzamos rcx (longitud del k-mer) posiciones en kmerList
    dec rbx                      ; Disminuimos el número de k-mers restantes por copiar
    jmp .copiar                  ; Repetimos el proceso para el siguiente k-mer

.fin:
    ret                          ; Terminamos y regresamos
;-------------------------------------------------------
; Ordenar k-mers usando bubble sort
ordenarKmers:
    push rbx              ; Guardamos todos los registros que vamos a usar
    push rdi
    push rsi
    push rcx
    push rdx
    push r8
    push r9
    push r12

    movzx r8, byte [kVal] ; r8 = longitud de cada k-mer
    mov rbx, [totalKmers] ; rbx = número total de k-mers

    cmp rbx, 1            ; Si hay 0 o 1 k-mer, ya está ordenado
    jle .fin
    dec rbx               ; rbx = totalKmers - 1 → número de pasadas necesarias

.outer:
    mov rcx, rbx          ; rcx = número de comparaciones en esta pasada
    xor r12, r12          ; r12 = índice base (offset) para el primer k-mer

.inner:
    mov rsi, r12
    add rsi, r8           ; rsi = índice del segundo k-mer a comparar (siguiente en lista)
    mov r9, 0             ; r9 = contador para comparar carácter por carácter del k-mer

.compare_loop:
    cmp r9, r8            ; ¿Ya se compararon todos los caracteres del k-mer?
    je .no_swap           ; Si sí, no se hace intercambio (son iguales)
    
    mov al, [kmerList + r12 + r9]  ; al = letra de primer k-mer
    mov dl, [kmerList + rsi + r9]  ; dl = letra de segundo k-mer
    cmp al, dl
    jb .no_swap           ; Si al < dl, ya están en orden → no swap
    ja .do_swap           ; Si al > dl, deben intercambiarse
    inc r9                ; Siguiente carácter
    jmp .compare_loop     ; Repetir comparación

.do_swap:
    xor r9, r9            ; r9 = 0 → usarlo para intercambiar cada letra del k-mer

.swap_loop:
    cmp r9, r8            ; ¿Ya se intercambiaron todos los caracteres?
    je .next              ; Si sí, saltar a siguiente par de k-mers
    mov al, [kmerList + r12 + r9] ; al = letra del primer k-mer
    mov dl, [kmerList + rsi + r9] ; dl = letra del segundo k-mer
    mov [kmerList + r12 + r9], dl ; intercambiar
    mov [kmerList + rsi + r9], al
    inc r9
    jmp .swap_loop        ; Repetir para todos los caracteres del k-mer

.no_swap:
.next:
    add r12, r8           ; Avanzar al siguiente k-mer
    dec rcx               ; Una comparación menos en esta pasada
    jnz .inner            ; Repetir bucle interno mientras queden comparaciones
    dec rbx               ; Una pasada menos por hacer
    jnz .outer            ; Repetir bucle externo

.fin:
    ; Restaurar los registros
    pop r12
    pop r9
    pop r8
    pop rdx
    pop rcx
    pop rsi
    pop rdi
    pop rbx
    ret
;-------------------------------------------------------
; Contar frecuencias de k-mers ordenados
contar_frecuencias:
    mov rbx, [totalKmers]         ; rbx = número total de k-mers
    cmp rbx, 0                    ; Si no hay ninguno, termina
    je .fin

    xor rsi, rsi                  ; rsi = índice actual dentro de kmerList
    movzx r8, byte [kVal]         ; r8 = longitud de cada k-mer
    mov r9, 1                     ; r9 = contador de repeticiones del k-mer actual

    ; Inicializar con el primer k-mer (se copia a tempKmer)
    call limpiar_tempKmer
    xor rcx, rcx
.copy_first:
    cmp rcx, r8                   ; ¿Ya copiamos todos los caracteres del primer k-mer?
    je .start_loop
    mov al, [kmerList + rsi + rcx] ; al = carácter del k-mer en la lista
    mov [tempKmer + rcx], al     ; Guardarlo en tempKmer
    inc rcx
    jmp .copy_first

.start_loop:
    add rsi, r8                   ; rsi apunta al siguiente k-mer
    dec rbx                       ; Ya procesamos uno
    jz .write_last                ; Si solo había uno, lo escribimos y terminamos

.loop:
    ; Comparar el tempKmer con el siguiente en la lista
    xor rcx, rcx
    mov r10, 1                    ; Asumimos que son iguales (bandera)

.compare_loop:
    cmp rcx, r8                   ; ¿Ya se compararon todos los caracteres?
    je .compare_done
    mov al, [tempKmer + rcx]      ; Carácter actual de tempKmer
    mov dl, [kmerList + rsi + rcx]; Carácter actual del k-mer nuevo
    cmp al, dl
    jne .different
    inc rcx
    jmp .compare_loop

.different:
    mov r10, 0                    ; Si un carácter difiere, ya no son iguales

.compare_done:
    cmp r10, 1                    ; ¿Son iguales?
    je .same_kmer                 ; Si sí, solo incrementar contador

    ; Si son diferentes:
    call escribir_kmer_freq       ; Escribir tempKmer con su frecuencia

    ; Copiar el nuevo k-mer a tempKmer
    call limpiar_tempKmer
    xor rcx, rcx
.copy_new:
    cmp rcx, r8
    je .reset_counter
    mov al, [kmerList + rsi + rcx]
    mov [tempKmer + rcx], al
    inc rcx
    jmp .copy_new

.reset_counter:
    mov r9, 1                     ; Reiniciar contador para nuevo k-mer
    jmp .continue

.same_kmer:
    inc r9                        ; Si es el mismo k-mer, incrementar su contador

.continue:
    add rsi, r8                   ; Pasamos al siguiente k-mer
    dec rbx                       ; Disminuimos el total restante
    jnz .loop                     ; Si quedan más, continuar

.write_last:
    ; Escribir el último grupo de k-mers iguales
    call escribir_kmer_freq

.fin:
    ret

;-------------------------------------------------------
; Escribir k-mer y su frecuencia en archivo
escribir_kmer_freq:
    push rax                    ; Guardamos registros usados
    push rbx
    push rcx
    push rdx
    push rsi
    push rdi
    
    movzx r8, byte [kVal]       ; r8 = longitud actual de los k-mers (de 4 a 10)

    ; Escribir el k-mer (contenido de tempKmer)
    mov rax, SYS_write          ; syscall write
    mov rdi, [outFileDesc]      ; descriptor del archivo de salida
    mov rsi, tempKmer           ; puntero al k-mer actual
    mov rdx, r8                 ; longitud del k-mer
    syscall                     ; escribe el k-mer en el archivo

    ; Escribir el separador ": "
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, separator          ; cadena ": "
    mov rdx, 2                  ; longitud de ": "
    syscall                     ; escribe el separador

    ; Escribir la frecuencia (guardada en r9)
    mov rax, r9                 ; rax = frecuencia del k-mer
    call print_decimal_to_file ; convertir y escribir número en el archivo

    ; Escribir salto de línea (crlf = \r\n)
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, crlf               ; \r\n
    mov rdx, 2
    syscall

    ; Restaurar los registros usados
    pop rdi
    pop rsi
    pop rdx
    pop rcx
    pop rbx
    pop rax
    ret                         ; regresar a la función que llamó
;-------------------------------------------------------
; Convertir número en rax a decimal y escribir en archivo
print_decimal_to_file:
    push rbx                     ; Guardar registros usados
    push rcx
    push rdx
    push rsi

    mov rdx, rax                 ; Guardamos el número original en rdx

    ; Limpiar 20 bytes del final del buffer (para asegurar espacio limpio)
    mov rcx, 20                 ; rcx será el contador
    mov rax, BUFF_SIZE         ; rax = 65536 (fin del buffer)
    sub rax, rcx               ; rax = inicio de los últimos 20 bytes
    lea rsi, [readBuffer + rax] ; rsi apunta a ese espacio en el buffer

.clear_loop:
    mov byte [rsi], 0          ; limpiamos byte por byte
    inc rsi
    loop .clear_loop           ; repetir 20 veces

    ; Empezamos desde el final del buffer para escribir el número al revés
    lea rbx, [readBuffer + BUFF_SIZE] ; rbx apunta justo al final del buffer
    mov rcx, 10                       ; divisor para decimal
    cmp rdx, 0                        ; ¿es cero el número?
    jne .convert
    dec rbx                           ; dejar espacio para un carácter
    mov byte [rbx], '0'              ; guardar el carácter '0'
    jmp .print                        ; saltar a impresión

.convert:
    xor rax, rax                     ; limpiar rax para la división
    mov rax, rdx                     ; rax = número a convertir
    xor rdx, rdx                     ; limpiar rdx

.div_loop:
    div rcx                          ; divide rax entre 10 → cociente en rax, residuo en rdx
    add dl, '0'                      ; convertir dígito a carácter ASCII
    dec rbx                          ; retroceder un byte en el buffer
    mov [rbx], dl                    ; guardar carácter
    xor rdx, rdx                     ; limpiar rdx para siguiente división
    cmp rax, 0
    jne .div_loop                    ; repetir hasta que rax (el número) sea 0

.print:
    mov rsi, rbx                     ; rsi = inicio del número en el buffer
    mov rdx, readBuffer + BUFF_SIZE ; calcular la longitud a escribir
    sub rdx, rbx                     ; longitud = fin - inicio
    mov rax, SYS_write
    mov rdi, [outFileDesc]          ; descriptor del archivo de salida
    syscall                         ; escribir el número convertido

    ; Restaurar registros
    pop rsi
    pop rdx
    pop rcx
    pop rbx
    ret                             ; regresar al que llamó

;-------------------------------------------
; Abrir archivo de entrada
openInputFile:
    mov rax, SYS_open          ; syscall para abrir archivo
    mov rdi, fileName          ; nombre del archivo a abrir (puntero a string)
    mov rsi, O_RDONLY          ; modo de solo lectura
    syscall                    ; ejecutar syscall (devuelve descriptor en rax)
    cmp rax, 0                 ; ¿se abrió correctamente?
    jl .error                  ; si rax < 0, hubo error
    mov [fileDesc], rax        ; guardar descriptor en variable fileDesc
    xor rax, rax               ; retornar 0 como éxito
    ret

.error:
    mov rdi, errMsgOpen        ; cargar mensaje de error
    call printStr              ; imprimir el mensaje
    mov rax, -1                ; retornar -1 como error
    ret


;-------------------------------------------
; Leer archivo FASTA
readFASTA:
    xor r15, r15               ; índice para llenar adnBuffer (contador de ADN)
    xor r12, r12               ; flag para saber si estamos dentro de una cabecera

.nextBlock:
    mov rax, SYS_read          ; syscall para leer archivo
    mov rdi, [fileDesc]        ; descriptor del archivo
    mov rsi, readBuffer        ; buffer donde leeremos los datos
    mov rdx, BUFF_SIZE         ; tamaño de lectura
    syscall                    ; ejecutar lectura
    cmp rax, 0                 ; ¿Fin de archivo?
    jle .done                  ; si <= 0, termina lectura
    mov r14, rax               ; r14 = bytes leídos
    xor r13, r13               ; r13 = índice dentro del bloque leído

.nextByte:
    cmp r13, r14               ; ¿ya recorrimos todo el bloque?
    je .nextBlock              ; si sí, leer siguiente bloque
    mov al, [readBuffer + r13] ; al = siguiente byte del bloque
    inc r13                    ; avanzar al siguiente byte

    test r12, r12              ; ¿estamos dentro de una cabecera?
    jnz .inHeader              ; si sí, seguir saltando hasta LF

    cmp al, '>'                ; ¿es inicio de cabecera?
    je .setHeader              ; marcar que estamos en cabecera
    cmp al, LF                 ; ¿es salto de línea?
    je .skip                   ; ignorar
    cmp al, 13                 ; ¿retorno de carro?
    je .skip                   ; ignorar
    cmp al, ';'                ; comentario en FASTA
    je .skip                   ; ignorar

    ; Validar letras válidas de ADN
    cmp al, 'A'
    je .store
    cmp al, 'C'
    je .store
    cmp al, 'G'
    je .store
    cmp al, 'T'
    je .store  

    ; Si es una letra minúscula, convertirla a mayúscula
    cmp al, 'a'
    jb .skip                   ; si menor que 'a', ignorar
    cmp al, 'z'
    ja .skip                   ; si mayor que 'z', ignorar
    sub al, 32                 ; convertir a mayúscula (ASCII trick)

.store:
    mov [adnBuffer + r15], al ; guardar la letra en adnBuffer
    inc r15                   ; incrementar índice de adnBuffer
    cmp r15, MAX_ADN_SIZE     ; ¿se alcanzó el límite?
    jae .done                 ; si sí, salir
.skip:
    jmp .nextByte             ; procesar siguiente byte

.setHeader:
    mov r12, 1                ; marcar que estamos dentro de cabecera
    jmp .skip                 ; ignorar este byte

.inHeader:
    cmp al, LF                ; ¿ya terminó la cabecera?
    jne .skip                 ; si no, seguir ignorando
    xor r12, r12              ; salir del modo cabecera
    jmp .skip

.done:
    mov [adnLength], r15      ; guardar longitud total de ADN válido
    xor rax, rax              ; retornar 0 como éxito
    ret
;-------------------------------------------
; Solicitar valor de k
getKValue:
    mov rdi, promptK          ; Carga en rdi el mensaje "Ingrese longitud de palabra..."
    call printStr             ; Llama a la función que imprime cadenas por consola

    ; Leer hasta 4 bytes para capturar números de 1-2 dígitos + enter
    mov rax, SYS_read         ; Preparar syscall para leer
    mov rdi, STDIN            ; Leer desde entrada estándar (teclado)
    mov rsi, readBuffer       ; Guardar la entrada en readBuffer
    mov rdx, 4                ; Leer máximo 4 bytes (ej. "10\n")
    syscall                   ; Ejecutar la lectura

    ; Verificar si el primer carácter es un dígito
    movzx rax, byte [readBuffer] ; Cargar el primer carácter como byte extendido a rax
    cmp rax, '0'              ; ¿es menor que '0'?
    jl .error                 ; Si sí, no es un dígito => error
    cmp rax, '9'              ; ¿es mayor que '9'?
    jg .error                 ; Si sí, no es un dígito => error

    ; Convertir primer dígito
    sub rax, '0'              ; Convierte ASCII a valor numérico (ej. '4' -> 4)
    mov rbx, rax              ; Guardar el primer dígito convertido en rbx

    ; Verificar si hay un segundo dígito
    movzx rax, byte [readBuffer + 1] ; Cargar segundo carácter (si existe)
    cmp rax, '0'
    jl .single_digit          ; Si no es dígito, usamos solo el primero
    cmp rax, '9'
    jg .single_digit          ; Si no es dígito, usamos solo el primero
    cmp rax, 10               ; ¿Es salto de línea?
    je .single_digit
    cmp rax, 13               ; ¿Es retorno de carro?
    je .single_digit

    ; Hay segundo dígito válido
    sub rax, '0'              ; Convertir segundo dígito de ASCII a número
    imul rbx, 10              ; Multiplicar el primero por 10 (para desplazar decimal)
    add rbx, rax              ; Sumar el segundo dígito (ej. 1 * 10 + 0 = 10)
    jmp .validate

.single_digit:
    ; Solo hay un dígito, rbx ya contiene el valor

.validate:
    ; Validar que esté en el rango 4-10
    cmp rbx, 4
    jl .error                 ; Si es menor a 4 => error
    cmp rbx, 10
    jg .error                 ; Si es mayor a 10 => error

    ; Guardar el valor válido
    mov [kVal], bl            ; Guardar el valor final de k en la variable kVal
    xor rax, rax              ; Retornar 0 como éxito
    ret

.error:
    mov rdi, errMsgK          ; Cargar mensaje de error
    call printStr             ; Imprimir mensaje de error
    mov rax, -1               ; Retornar -1 como error
    ret

;-------------------------------------------
; Solo procesa los k-mers, no los escribe desordenados
generateWords:
    ; Solo retornamos éxito, el procesamiento real se hace después
    xor rax, rax       ; Coloca 0 en rax para indicar éxito (sin error)
    ret                ; Regresa de la función

;-------------------------------------------------------
; Limpiar el buffer tempKmer antes de cada uso
limpiar_tempKmer:
    push rax           ; Guarda rax en la pila para no perder su valor
    push rcx           ; Guarda rcx en la pila para no perder su valor
    xor rcx, rcx       ; rcx = 0 (contador para el bucle)
    mov rax, 16        ; rax = 16 (longitud del buffer a limpiar)
.limpiar_loop:
    mov byte [tempKmer + rcx], 0 ; Escribe 0 en la posición actual de tempKmer
    inc rcx            ; Incrementa el contador
    dec rax            ; Decrementa la cantidad restante a limpiar
    jnz .limpiar_loop  ; Si no es cero, repite el bucle
    pop rcx            ; Restaura el valor original de rcx
    pop rax            ; Restaura el valor original de rax
    ret                ; Regresa de la función
;-------------------------------------------------------
; Final del código