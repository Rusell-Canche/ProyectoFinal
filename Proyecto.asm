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
crlf db 13, 10, 0

section .bss
readBuffer resb BUFF_SIZE

adnChunk resb MAX_ADN_CHUNK + 16     ; chunk con margen
overlapSize resq 1

fileDesc resq 1
outFileDesc resq 1
kVal resb 1
adnLength resq 1
kmerList resb MAX_ADN_SIZE
kmerCounts resq MAX_ADN_SIZE
tempKmer resb 16
contador resq 1
totalKmers resq 1

section .text
global _start

_start:
    mov rdi, header
    call printStr

    call openInputFile
    cmp rax, 0
    jl _exitError

    call getKValue
    cmp rax, 0
    jl _exitError

    ; Leer archivo FASTA
    call readFastaFile
    cmp rax, 0
    jl _exitError

    ; Cerrar archivo de entrada
    mov rax, SYS_close
    mov rdi, [fileDesc]
    syscall

    call generateWords
    cmp rax, 0
    jl _exitError

    ; Crear archivo de salida (truncar si existe)
    mov rax, SYS_creat
    mov rdi, outFileName
    mov rsi, S_IRUSR | S_IWUSR
    syscall
    cmp rax, 0
    jl _exitError
    mov [outFileDesc], rax

    call ordenarYGuardarKmers
    call ordenarKmers

    mov rdi, ordenadoMsg
    call printStrToFile

    call contar_frecuencias

    ; Cerrar archivo de salida
    mov rax, SYS_close
    mov rdi, [outFileDesc]
    syscall

    mov rdi, successMsg
    call printStr

    mov rax, SYS_exit
    mov rdi, EXIT_SUCCESS
    syscall
;-------------------------------------------

_exitError:
    mov rax, SYS_exit
    mov rdi, 1
    syscall

;-------------------------------------------
; Función para imprimir cadenas por consola
printStr:
    push rbx
    mov rbx, rdi
    xor rdx, rdx
.cont:
    cmp byte [rbx + rdx], 0
    je .done
    inc rdx
    jmp .cont
.done:
    mov rax, SYS_write
    mov rdi, STDOUT
    mov rsi, rbx
    syscall
    pop rbx
    ret

;-------------------------------------------
; Función para imprimir cadenas en archivo
printStrToFile:
    push rbx
    mov rbx, rdi
    xor rdx, rdx
.lenloop:
    cmp byte [rbx + rdx], 0
    je .write
    inc rdx
    jmp .lenloop
.write:
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, rbx
    syscall
    pop rbx
    ret

;-------------------------------------------
; Leer archivo FASTA/FNA
readFastaFile:
    push rbx
    push rcx
    push rdx
    push rsi
    push rdi

    xor rbx, rbx                ; índice en adnChunk
    
.read_loop:
    ; Leer chunk del archivo
    mov rax, SYS_read
    mov rdi, [fileDesc]
    mov rsi, readBuffer
    mov rdx, BUFF_SIZE
    syscall
    
    cmp rax, 0
    jle .done_reading           ; EOF o error
    
    mov rcx, rax                ; bytes leídos
    xor rsi, rsi                ; índice en readBuffer
    
.process_chunk:
    cmp rsi, rcx
    jge .read_loop
    
    mov al, [readBuffer + rsi]
    
    ; Saltar líneas de cabecera (que empiezan con >)
    cmp al, '>'
    je .skip_header_line
    
    ; Saltar caracteres de control
    cmp al, 10                  ; LF
    je .next_char
    cmp al, 13                  ; CR
    je .next_char
    cmp al, 32                  ; espacio
    je .next_char
    cmp al, 9                   ; tab
    je .next_char
    
    ; Validar que sea nucleótido válido (A, T, G, C)
    cmp al, 'A'
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
    jmp .next_char              ; Saltar caracteres no válidos
    
.convert_to_upper:
    sub al, 32                  ; Convertir a mayúscula
    
.valid_nucleotide:
    ; Verificar que no excedamos el buffer
    cmp rbx, MAX_ADN_CHUNK - 1
    jge .done_reading
    
    mov [adnChunk + rbx], al
    inc rbx
    jmp .next_char
    
.skip_header_line:
    ; Saltar hasta el final de la línea
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
    mov [adnLength], rbx        ; Guardar longitud total
    
    ; Verificar que tengamos datos
    cmp rbx, 0
    je .error
    
    xor rax, rax                ; Éxito
    jmp .exit
    
.error:
    mov rdi, errMsgRead
    call printStr
    mov rax, -1
    
.exit:
    pop rdi
    pop rsi
    pop rdx
    pop rcx
    pop rbx
    ret

;-------------------------------------------------------
; Copiar k-mers en kmerList para ordenarlos
ordenarYGuardarKmers:
    xor rsi, rsi                  ; índice dentro del chunk
    xor rdi, rdi                  ; índice dentro de kmerList
    movzx rcx, byte [kVal]       ; longitud del k-mer
    mov rbx, [adnLength]         ; bytes válidos en adnChunk

    cmp rbx, rcx
    jb .fin

    sub rbx, rcx                 ; número de k-mers posibles
    inc rbx
    mov [totalKmers], rbx        ; guardar total de k-mers

.copiar:
    cmp rbx, 0
    je .fin

    mov rdx, 0
.loop_kmer:
    cmp rdx, rcx
    je .next
    mov al, [adnChunk + rsi + rdx]
    mov [kmerList + rdi + rdx], al
    inc rdx
    jmp .loop_kmer

.next:
    inc rsi
    add rdi, rcx
    dec rbx
    jmp .copiar

.fin:
    ret

;-------------------------------------------------------
; Ordenar k-mers usando bubble sort
ordenarKmers:
    push rbx
    push rdi
    push rsi
    push rcx
    push rdx
    push r8
    push r9
    push r10
    push r11
    push r12

    movzx r8, byte [kVal]      ; tamaño de k-mer
    mov rbx, [totalKmers]
    cmp rbx, 1
    jle .fin
    dec rbx                    ; n - 1 iteraciones

.outer:
    mov rcx, rbx
    xor r12, r12               ; offset inicial
    mov r11b, 0                ; flag de intercambio

.inner:
    mov rsi, r12
    add rsi, r8                ; siguiente k-mer
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
    ; Si hubo intercambio, marcarlo
    mov r11b, 1

.no_swap:
.next:
    add r12, r8
    dec rcx
    jnz .inner

    ; Si no hubo intercambios, está ordenado
    cmp r11b, 0
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
; Contar frecuencias de k-mers ordenados
contar_frecuencias:
    mov rbx, [totalKmers]
    cmp rbx, 0
    je .fin

    xor rsi, rsi               ; índice actual en kmerList
    movzx r8, byte [kVal]      ; longitud del k-mer
    mov r9, 1                  ; contador de repeticiones del k-mer actual

    ; Copiar primer k-mer a tempKmer
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
    jz .write_last             ; solo había un elemento

.loop:
    ; Comparar tempKmer con kmer actual
    xor rcx, rcx
    mov r10, 1                 ; flag = iguales

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

    ; Son diferentes: escribir anterior y copiar nuevo
    call escribir_kmer_freq

    ; Copiar nuevo a tempKmer
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
; Escribir k-mer y su frecuencia en archivo
escribir_kmer_freq:
    push rax
    push rbx
    push rcx
    push rdx
    push rsi
    push rdi
    
    movzx r8, byte [kVal]

    ; Escribir el k-mer
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, tempKmer
    mov rdx, r8
    syscall

    ; Escribir el separador
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, separator
    mov rdx, 2
    syscall

    ; Escribir la frecuencia
    mov rax, r9
    call print_decimal_to_file

    ; Escribir salto de línea
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
; Convertir número en rax a decimal y escribir en archivo
print_decimal_to_file:
    push rbx
    push rcx
    push rdx
    push rsi

    mov rdx, rax                ; Guardamos el número original

    ; Limpiar 20 bytes del final del buffer
    mov rcx, 20
    mov rax, BUFF_SIZE
    sub rax, rcx
    lea rsi, [readBuffer + rax]

.clear_loop:
    mov byte [rsi], 0
    inc rsi
    loop .clear_loop

    ; Empezamos desde el final del buffer
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
; Abrir archivo de entrada
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
; Solicitar valor de k
getKValue:
    mov rdi, promptK
    call printStr

    mov rax, SYS_read
    mov rdi, STDIN
    mov rsi, readBuffer
    mov rdx, 4          ; leer hasta 4 caracteres incluyendo '\n'
    syscall

    xor rcx, rcx        ; índice
    xor rbx, rbx        ; acumulador para el número

.parse_loop:
    mov al, [readBuffer + rcx]
    cmp al, 10          ; si es salto de línea (\n)
    je .done_parse
    cmp al, 13          ; si es retorno de carro (por si acaso)
    je .done_parse
    cmp al, 0
    je .done_parse
    sub al, '0'
    cmp al, 9
    ja .error           ; si no es dígito

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
; Generar k-mers (ahora implementada correctamente)
generateWords:
    ; Esta función ahora solo retorna éxito ya que el procesamiento
    ; real se hace en ordenarYGuardarKmers
    xor rax, rax
    ret

;-------------------------------------------------------
; Limpiar el buffer tempKmer antes de cada uso
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