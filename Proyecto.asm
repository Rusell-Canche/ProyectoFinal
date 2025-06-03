; Proyecto.asm
; Programa para analizar palabras de ADN en archivos FASTA/FNA
; Características:
;   - Lee archivos FASTA/FNA, omite líneas de cabecera
;   - Genera palabras de longitud k (4-10)
;   - Guarda SOLO palabras ordenadas con repeticiones en "palabras.txt"

section .data
LF equ 10
NULL equ 0
EXIT_SUCCESS equ 0
STDIN equ 0
STDOUT equ 1
STDERR equ 2
SYS_read equ 0
SYS_write equ 1
SYS_open equ 2
SYS_close equ 3
SYS_exit equ 60
SYS_creat equ 85
O_RDONLY equ 000000q
O_WRONLY equ 000001q
O_CREAT equ 0x40
O_TRUNC equ 0x200
O_APPEND equ 0x400

S_IRUSR equ 00400q
S_IWUSR equ 00200q
BUFF_SIZE equ 65536     ; 64KB por lectura

MAX_ADN_CHUNK  equ 524288     ; 512 KB por bloque
MAX_ADN_SIZE   equ 524288     ; Definir tamaño máximo para kmerList y kmerCounts

separator db ": ", 0

header db "Analizador de palabras de ADN", LF, LF, NULL
fileName db "Prueba.txt", NULL
outFileName db "palabras.txt", NULL
promptK db "Ingrese longitud de palabra (4-10): ", NULL
espacio db " ", 0
ordenadoMsg db "Palabras Ordenadas con Repeticiones:", 13, 10, 0
repeticionesMsg db "Repeticiones:", 13, 10, 0

errMsgOpen db "Error abriendo archivo.", LF, NULL
errMsgRead db "Error leyendo archivo.", LF, NULL
errMsgWrite db "Error escribiendo archivo.", LF, NULL
errMsgK db "Error: k debe ser entre 4 y 10.", LF, NULL
successMsg db "Palabras guardadas en palabras.txt", LF, NULL
crlf db 13, 10, 0   ; Carácter de salto de línea (CRLF)

section .bss
readBuffer resb BUFF_SIZE         ; Buffer para leer del archivo
adnChunk resb MAX_ADN_CHUNK + 16  ; Buffer para almacenar la secuencia de ADN leída
overlapSize resq 1                ; No usado, reservado para posibles solapamientos
fileDesc resq 1                   ; Descriptor de archivo de entrada
outFileDesc resq 1                ; Descriptor de archivo de salida
kVal resb 1                       ; Valor de k (longitud de las palabras)
adnLength resq 1                  ; Longitud de la secuencia de ADN leída
kmerList resb MAX_ADN_SIZE        ; Lista de k-mers extraídos
kmerCounts resq MAX_ADN_SIZE      ; No usado, reservado para contar repeticiones
tempKmer resb 16                  ; Buffer temporal para comparar k-mers
contador resq 1                   ; No usado, reservado para contar
totalKmers resq 1                 ; Total de k-mers generados

section .text
global _start

;-------------------------------------------
; Punto de entrada principal del programa
_start:
    mov rdi, header
    call printStr                ; Imprime el encabezado

    call openInputFile           ; Abre el archivo de entrada
    cmp rax, 0                   ; Verifica si hubo error al abrir
    jl _exitError                ; Si hay error, termina

    call getKValue               ; Solicita el valor de k al usuario
    cmp rax, 0                   ; Verifica si hubo error al obtener k
    jl _exitError

    call readFastaFile           ; Lee el archivo FASTA y almacena la secuencia de ADN
    cmp rax, 0                   ; Verifica si hubo error al leer
    jl _exitError

    ; Cierra el archivo de entrada
    mov rax, SYS_close           ; Cierra el archivo de entrada
    mov rdi, [fileDesc]          ; descriptor del archivo
    syscall                      ;Sirve para cerrar el archivo

    call generateWords           ; (Actualmente solo retorna éxito)

    cmp rax, 0                   ; Verifica si hubo error al generar palabras
    jl _exitError                ; Si hay error, termina

    ; Crea el archivo de salida (trunca si existe)
    mov rax, SYS_creat           ; Crea o trunca el archivo de salida
    mov rdi, outFileName         ; Nombre del archivo de salida
    mov rsi, S_IRUSR | S_IWUSR   ; Permisos de lectura y escritura para el usuario
    syscall                      ; Llama al sistema para crear el archivo
    cmp rax, 0                   ; Verifica si hubo error al crear el archivo
    jl _exitError                ; Si hay error, termina
    mov [outFileDesc], rax       ; Guarda el descriptor del archivo de salida

    call ordenarYGuardarKmers    ; Extrae y guarda los k-mers en kmerList
    call ordenarKmers            ; Ordena los k-mers alfabéticamente

    mov rdi, ordenadoMsg         ; Mensaje de k-mers ordenados
    call printStrToFile          ; Escribe mensaje de encabezado en el archivo de salida

    call contar_frecuencias      ; Cuenta y escribe las frecuencias de cada k-mer

    ; Cierra el archivo de salida
    mov rax, SYS_close
    mov rdi, [outFileDesc]        ; descriptor del archivo de salida
    syscall                       ; Cierra el archivo de salida

    mov rdi, successMsg           ; Mensaje de éxito
    call printStr                 ; Mensaje de éxito

    mov rax, SYS_exit             ; Termina el programa con éxito
    mov rdi, EXIT_SUCCESS         ;sirve para indicar que el programa terminó correctamente
    syscall                       ;sirve para terminar el programa

;-------------------------------------------
; Termina el programa con error
_exitError:
    mov rax, SYS_exit              ; Termina el programa con error
    mov rdi, 1                     ; Código de error 1
    syscall                        ; Cierra el archivo de salida si está abierto

;-------------------------------------------
; Imprime una cadena por consola (rdi = puntero a la cadena)
printStr:
    push rbx                       ; sirve para guardar el registro rbx
    mov rbx, rdi                   ; rbx apunta a la cadena a imprimir
    xor rdx, rdx                   ;comparador de longitud de cadena
.cont:
    cmp byte [rbx + rdx], 0        ; verifica si el final de la cadena es nulo 
    je .done                       ; si es nulo, termina
    inc rdx                        ; incrementa el contador de longitud                      
    jmp .cont                      ; vuelve al inicio del bucle
.done:
    mov rax, SYS_write             ; Llama al sistema para escribir
    mov rdi, STDOUT                ; Descriptor de archivo estándar de salida
    mov rsi, rbx                   ; Puntero a la cadena
    syscall                        ; Llama al sistema para escribir la cadena
    pop rbx                        ; Recupera el registro rbx
    ret                            ;sirve para retornar al punto de llamada

;-------------------------------------------
; Imprime una cadena en el archivo de salida (rdi = puntero a la cadena)
printStrToFile:
    push rbx                        ; Guarda el registro rbx
    mov rbx, rdi                    ; rbx apunta a la cadena a imprimir
    xor rdx, rdx                    ; Contador de longitud de cadena
.lenloop:
    cmp byte [rbx + rdx], 0         ; Verifica si el final de la cadena es nulo
    je .write                       ; Si es nulo, escribe
    inc rdx                         ; Incrementa el contador de longitud
    jmp .lenloop
.write:
    mov rax, SYS_write              ; Llama al sistema para escribir
    mov rdi, [outFileDesc]          ; Descriptor de archivo de salida
    mov rsi, rbx                    ; Puntero a la cadena
    syscall                         ; Llama al sistema para escribir la cadena
    pop rbx                         ; Recupera el registro rbx para no perder su valor
    ret

;-------------------------------------------
; Lee el archivo FASTA/FNA y almacena la secuencia de ADN en adnChunk
readFastaFile:
    push rbx                         ; Guarda el registro rbx que es usado como índice de lectura
    push rcx                         ; Guarda el registro rcx que es usado para contar bytes leídos
    push rdx                         ; Guarda el registro rdx que es usado para manejar el buffer de lectura
    push rsi                         ; Guarda el registro rsi que es usado como índice en readBuffer
    push rdi                         ; Guarda el registro rdi que es usado para manejar el descriptor de archivo

    xor rbx, rbx                     ; Índice en adnChunk que se va a llenar
    
.read_loop:
    ; Lee un chunk del archivo
    mov rax, SYS_read                 ; Llama al sistema para leer
    mov rdi, [fileDesc]               ; Descriptor del archivo de entrada
    mov rsi, readBuffer               ; Puntero al buffer de lectura
    mov rdx, BUFF_SIZE                ; Tamaño del buffer
    syscall
    
    cmp rax, 0
    jle .done_reading           ; EOF o error
    
    mov rcx, rax                ; Bytes leídos
    xor rsi, rsi                ; Índice en readBuffer
    
.process_chunk:
    cmp rsi, rcx                ; Verifica si se han procesado todos los bytes leídos
    jge .read_loop              ; Si sí, lee el siguiente chunk
    
    mov al, [readBuffer + rsi]  ; Carga el siguiente byte del buffer
    
    ; Salta líneas de cabecera (empiezan con '>')
    cmp al, '>'
    je .skip_header_line
    
    ; Salta caracteres de control (LF, CR, espacio, tab)
    cmp al, 10                  ; LF
    je .next_char
    cmp al, 13                  ; CR
    je .next_char
    cmp al, 32                  ; espacio
    je .next_char
    cmp al, 9                   ; tab
    je .next_char
    
    ; Valida que sea nucleótido válido (A, T, G, C, a, t, g, c)
    cmp al, 'A'                  ; Verifica si es un nucleótido válido
    je .valid_nucleotide
    cmp al, 'T' 
    je .valid_nucleotide
    cmp al, 'G'
    je .valid_nucleotide
    cmp al, 'C'
    je .valid_nucleotide
    cmp al, 'a'
    je .convert_to_upper
    cmp al, 't'
    je .convert_to_upper
    cmp al, 'g'
    je .convert_to_upper
    cmp al, 'c'
    je .convert_to_upper
    jmp .next_char              ; Salta caracteres no válidos
    
.convert_to_upper:
    sub al, 32                  ; Convierte a mayúscula
    
.valid_nucleotide:
    ; Verifica que no se exceda el buffer
    cmp rbx, MAX_ADN_CHUNK - 1
    jge .done_reading
    
    mov [adnChunk + rbx], al    ; Guarda el nucleótido
    inc rbx
    jmp .next_char
    
.skip_header_line:
    ; Salta hasta el final de la línea de cabecera
.skip_loop:
    inc rsi
    cmp rsi, rcx
    jge .read_loop
    mov al, [readBuffer + rsi]
    cmp al, 10                  ; LF
    jne .skip_loop
    
.next_char:
    inc rsi
    jmp .process_chunk
    
.done_reading:
    mov [adnLength], rbx        ; Guarda la longitud total de la secuencia
    
    ; Verifica que haya datos
    cmp rbx, 0
    je .error
    
    xor rax, rax                ; Éxito
    jmp .exit
    
.error:
    mov rdi, errMsgRead
    call printStr
    mov rax, -1
    
.exit:
    pop rdi                      ;sirve para recuperar el registro rdi
    pop rsi
    pop rdx
    pop rcx
    pop rbx
    ret

;-------------------------------------------------------
; Extrae todos los k-mers posibles de la secuencia y los guarda en kmerList
ordenarYGuardarKmers:
    xor rsi, rsi                  ; Índice dentro del chunk
    xor rdi, rdi                  ; Índice dentro de kmerList
    movzx rcx, byte [kVal]        ; Longitud del k-mer
    mov rbx, [adnLength]          ; Bytes válidos en adnChunk

    cmp rbx, rcx
    jb .fin                       ; Si la secuencia es menor que k, termina

    sub rbx, rcx                  ; Número de k-mers posibles
    inc rbx
    mov [totalKmers], rbx         ; Guarda el total de k-mers

.copiar:
    cmp rbx, 0
    je .fin

    mov rdx, 0
.loop_kmer:
    cmp rdx, rcx                  ; Verifica si se han copiado todos los caracteres del k-mer
    je .next
    mov al, [adnChunk + rsi + rdx] ; Carga el siguiente nucleótido
    mov [kmerList + rdi + rdx], al ; Guarda el nucleótido en kmerList
    inc rdx                       ; Incrementa el índice del k-mer
    jmp .loop_kmer

.next:
    inc rsi                       ; Avanza al siguiente nucleótido en adnChunk
    add rdi, rcx                  ; Avanza al siguiente k-mer en kmerList
    dec rbx                       ; Decrementa el contador de k-mers restantes
    jmp .copiar

.fin:
    ret

;-------------------------------------------------------
; Ordena los k-mers en kmerList usando bubble sort
ordenarKmers:
    push rbx                       ; Guarda el registro rbx que es usado como contador de k-mers
    push rdi                       ; Guarda el registro rdi que es usado como índice de kmerList
    push rsi                       ; Guarda el registro rsi que es usado como índice de comparación
    push rcx                       ; Guarda el registro rcx que es usado como contador de iteraciones
    push r8                        ; Guarda el registro r8 que es usado para el tamaño del k-mer
    push r9                        ; Guarda el registro r9 que es usado como índice de comparación
    push r10                       ; Guarda el registro r10 que es usado como flag de intercambio
    push r11                       ; Guarda el registro r11 que es usado como flag de intercambio
    push r12                       ; Guarda el registro r12 que es usado como offset inicial

    movzx r8, byte [kVal]      ; Tamaño de k-mer
    mov rbx, [totalKmers]
    cmp rbx, 1
    jle .fin
    dec rbx                    ; n - 1 iteraciones

.outer:
    mov rcx, rbx
    xor r12, r12               ; Offset inicial
    mov r11b, 0                ; Flag de intercambio

.inner:
    mov rsi, r12
    add rsi, r8                ; Siguiente k-mer
    mov r9, 0

.compare_loop:
    cmp r9, r8
    je .no_swap
    mov al, [kmerList + r12 + r9]
    mov dl, [kmerList + rsi + r9]
    cmp al, dl
    jb .no_swap
    ja .do_swap
    inc r9
    jmp .compare_loop

.do_swap:
    xor r9, r9
.swap_loop:
    cmp r9, r8
    je .swap_done
    mov al, [kmerList + r12 + r9]
    mov dl, [kmerList + rsi + r9]
    mov [kmerList + r12 + r9], dl
    mov [kmerList + rsi + r9], al
    inc r9
    jmp .swap_loop

.swap_done:
    mov r11b, 1                ; Marca que hubo intercambio

.no_swap:
.next:
    add r12, r8
    dec rcx
    jnz .inner

    cmp r11b, 0                ; Si no hubo intercambios, ya está ordenado
    je .fin
    dec rbx
    jnz .outer

.fin:
    pop r12
    pop r11
    pop r10
    pop r9
    pop r8
    pop rdx
    pop rcx
    pop rsi
    pop rdi
    pop rbx
    ret

;-------------------------------------------------------
; Cuenta las frecuencias de los k-mers ordenados y los escribe en el archivo
contar_frecuencias:
    mov rbx, [totalKmers]
    cmp rbx, 0
    je .fin

    xor rsi, rsi               ; Índice actual en kmerList
    movzx r8, byte [kVal]      ; Longitud del k-mer
    mov r9, 1                  ; Contador de repeticiones del k-mer actual

    ; Copia el primer k-mer a tempKmer
    call limpiar_tempKmer
    xor rcx, rcx
.copy_first:
    cmp rcx, r8
    je .start_loop
    mov al, [kmerList + rsi + rcx]
    mov [tempKmer + rcx], al
    inc rcx
    jmp .copy_first

.start_loop:
    add rsi, r8
    dec rbx
    jz .write_last             ; Solo había un elemento

.loop:
    ; Compara tempKmer con el k-mer actual
    xor rcx, rcx
    mov r10, 1                 ; Flag = iguales

.compare_loop:
    cmp rcx, r8
    je .compare_done
    mov al, [tempKmer + rcx]
    mov dl, [kmerList + rsi + rcx]
    cmp al, dl
    jne .not_equal
    inc rcx
    jmp .compare_loop

.not_equal:
    mov r10, 0

.compare_done:
    cmp r10, 1
    je .same_kmer

    ; Son diferentes: escribe el anterior y copia el nuevo
    call escribir_kmer_freq

    ; Copia el nuevo a tempKmer
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
    mov r9, 1
    jmp .continue

.same_kmer:
    inc r9

.continue:
    add rsi, r8
    dec rbx
    jnz .loop

.write_last:
    call escribir_kmer_freq

.fin:
    ret

;-------------------------------------------------------
; Escribe el k-mer y su frecuencia en el archivo de salida
escribir_kmer_freq:
    push rax
    push rbx
    push rcx
    push rdx
    push rsi
    push rdi
    
    movzx r8, byte [kVal]

    ; Escribe el k-mer
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, tempKmer
    mov rdx, r8
    syscall

    ; Escribe el separador ": "
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, separator
    mov rdx, 2
    syscall

    ; Escribe la frecuencia (número)
    mov rax, r9
    call print_decimal_to_file

    ; Escribe salto de línea
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, crlf
    mov rdx, 2
    syscall

    pop rdi
    pop rsi
    pop rdx
    pop rcx
    pop rbx
    pop rax
    ret

;-------------------------------------------------------
; Convierte el número en rax a decimal y lo escribe en el archivo de salida
print_decimal_to_file:
    push rbx
    push rcx
    push rdx
    push rsi

    mov rdx, rax                ; Guarda el número original

    ; Limpia 20 bytes del final del buffer
    mov rcx, 20
    mov rax, BUFF_SIZE
    sub rax, rcx
    lea rsi, [readBuffer + rax]

.clear_loop:
    mov byte [rsi], 0
    inc rsi
    loop .clear_loop

    ; Empieza desde el final del buffer
    lea rbx, [readBuffer + BUFF_SIZE]
    mov rcx, 10
    cmp rdx, 0
    jne .convert
    dec rbx
    mov byte [rbx], '0'
    jmp .print

.convert:
    xor rax, rax
    mov rax, rdx
    xor rdx, rdx
.div_loop:
    div rcx
    add dl, '0'
    dec rbx
    mov [rbx], dl
    xor rdx, rdx
    cmp rax, 0
    jne .div_loop

.print:
    mov rsi, rbx
    mov rdx, readBuffer + BUFF_SIZE
    sub rdx, rbx
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    syscall

    pop rsi
    pop rdx
    pop rcx
    pop rbx
    ret

;-------------------------------------------
; Abre el archivo de entrada y guarda el descriptor en fileDesc
openInputFile:
    mov rax, SYS_open
    mov rdi, fileName
    mov rsi, O_RDONLY
    syscall
    cmp rax, 0
    jl .error
    mov [fileDesc], rax
    xor rax, rax
    ret
.error:
    mov rdi, errMsgOpen
    call printStr
    mov rax, -1
    ret

;-------------------------------------------
; Solicita el valor de k al usuario y lo valida
getKValue:
    mov rdi, promptK
    call printStr

    mov rax, SYS_read
    mov rdi, STDIN
    mov rsi, readBuffer
    mov rdx, 4          ; Lee hasta 4 caracteres incluyendo '\n'
    syscall

    xor rcx, rcx        ; Índice
    xor rbx, rbx        ; Acumulador para el número

.parse_loop:
    mov al, [readBuffer + rcx]
    cmp al, 10          ; Si es salto de línea (\n)
    je .done_parse
    cmp al, 13          ; Si es retorno de carro
    je .done_parse
    cmp al, 0
    je .done_parse
    sub al, '0'
    cmp al, 9
    ja .error           ; Si no es dígito

    imul rbx, rbx, 10
    add rbx, rax

    inc rcx
    cmp rcx, 4
    jl .parse_loop

.done_parse:
    cmp rbx, 4
    jl .error
    cmp rbx, 10
    jg .error
    mov [kVal], bl
    xor rax, rax
    ret

.error:
    mov rdi, errMsgK
    call printStr
    mov rax, -1
    ret

;-------------------------------------------
; Función dummy: solo retorna éxito (el procesamiento real está en ordenarYGuardarKmers)
generateWords:
    xor rax, rax
    ret

;-------------------------------------------------------
; Limpia el buffer tempKmer antes de cada uso
limpiar_tempKmer:
    push rax
    push rcx
    xor rcx, rcx
    mov rax, 16        ; Tamaño máximo del buffer
.limpiar_loop:
    mov byte [tempKmer + rcx], 0
    inc rcx
    dec rax
    jnz .limpiar_loop
    pop rcx
    pop rax
    ret