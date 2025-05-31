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
section .data
    ; Constantes
    EXIT_SUCCESS    equ 0
    SYS_exit        equ 60
    SYS_read        equ 0
    SYS_write       equ 1
    SYS_open        equ 2
    SYS_close       equ 3
    SYS_creat       equ 85
    O_RDONLY        equ 0
    LF              equ 10
    NULL            equ 0
    STDOUT          equ 1
    STDIN           equ 0
    MAX_KMER_LEN    equ 10
    MAX_PATH_LEN    equ 256
    MAX_FILE_SIZE   equ 1000000
    MAX_KMERS       equ 1048576   ; 4^10 = 1,048,576

    ; Mensajes
    msg_k_prompt        db "Enter k (4-10): ", NULL
    msg_filename        db "Enter filename: ", NULL
    msg_invalid_k       db "Invalid k. Must be 4-10.", LF, NULL
    msg_file_error      db "Error opening file.", LF, NULL
    msg_plasmid         db "Plasmid detected. Sequence skipped.", LF, NULL
    msg_invalid_char    db "Invalid character found: ", NULL
    msg_clean_chars     db "Valid characters processed: ", NULL
    msg_total_kmers     db "Total k-mers generated: ", NULL
    msg_search_prompt   db "Search k-mer (or 'exit'): ", NULL
    msg_found           db "Frequency: ", NULL
    msg_not_found       db "K-mer not found.", LF, NULL
    msg_exit            db "Exiting...", LF, NULL
    newline             db LF, NULL
    tab                 db "    ", NULL

    ; Variables
    k_value         dq 0
    filename        times MAX_PATH_LEN db NULL
    file_desc       dq 0
    file_size       dq 0
    plasmid_flag    db 0
    invalid_count   dq 0
    clean_count     dq 0
    total_kmers     dq 0
    exit_cmd        db "exit", NULL

section .bss
    file_buffer     resb MAX_FILE_SIZE
    clean_dna       resb MAX_FILE_SIZE
    clean_len       resq 1
    freq_table      resd MAX_KMERS
    kmer_array      resb (16 * MAX_KMERS)  ; 12-byte string + 4-byte frequency
    num_non_zero    resq 1
    input_buffer    resb 256
    char_buffer     resb 2

section .text
    global _start

_start:
    ; Pedir valor de k
    mov rdi, msg_k_prompt
    call print_string
    mov rdi, input_buffer
    mov rsi, 256
    call read_input
    call parse_int
    cmp rax, 4
    jl .invalid_k
    cmp rax, 10
    jg .invalid_k
    mov [k_value], rax

    ; Pedir nombre de archivo
    mov rdi, msg_filename
    call print_string
    mov rdi, filename
    mov rsi, MAX_PATH_LEN
    call read_input

    ; Procesar archivo
    call open_file
    cmp rax, 0
    jl .file_error
    mov [file_desc], rax

    call read_file
    call close_file
    call process_dna
    call generate_kmers
    call build_kmer_array
    call quicksort
    jmp .search_loop

.invalid_k:
    mov rdi, msg_invalid_k
    call print_string
    jmp .exit

.file_error:
    mov rdi, msg_file_error
    call print_string
    jmp .exit

.search_loop:
    mov rdi, msg_search_prompt
    call print_string
    mov rdi, input_buffer
    mov rsi, 256
    call read_input

    ; Verificar si es "exit"
    mov rdi, input_buffer
    mov rsi, exit_cmd
    call strcmp
    cmp rax, 0
    je .exit

    ; Buscar k-mer
    mov rdi, input_buffer
    mov rsi, [k_value]
    call search_kmer
    cmp rax, -1
    je .not_found
    
    mov rdi, msg_found
    call print_string
    mov rdi, rax
    call print_uint
    mov rdi, newline
    call print_string
    jmp .search_loop

.not_found:
    mov rdi, msg_not_found
    call print_string
    jmp .search_loop

.exit:
    mov rdi, msg_exit
    call print_string
    mov rax, SYS_exit
    mov rdi, EXIT_SUCCESS
    syscall

; ===== Funciones principales =====

; Abrir archivo para lectura
open_file:
    mov rax, SYS_open
    mov rdi, filename
    mov rsi, O_RDONLY
    syscall
    ret

; Leer archivo completo
read_file:
    mov rax, SYS_read
    mov rdi, [file_desc]
    mov rsi, file_buffer
    mov rdx, MAX_FILE_SIZE
    syscall
    mov [file_size], rax
    ret

; Cerrar archivo
close_file:
    mov rax, SYS_close
    mov rdi, [file_desc]
    syscall
    ret

; Procesar secuencia de ADN
process_dna:
    mov rsi, file_buffer        ; Puntero al buffer del archivo
    mov rdi, clean_dna          ; Puntero al buffer limpio
    mov rcx, [file_size]        ; Tamaño del archivo
    xor rbx, rbx                ; Contador de caracteres limpios
    xor r8, r8                  ; Flag para plásmido (0 = no plásmido)

.process_loop:
    cmp rcx, 0
    je .process_done
    mov al, [rsi]

    ; Verificar si es inicio de cabecera
    cmp al, '>'
    jne .check_newline
    mov r8, 1                   ; Activar flag de plásmido
    jmp .skip_char

.check_newline:
    cmp al, LF
    jne .check_valid_char
    mov r8, 0                   ; Desactivar flag de plásmido al final de línea
    jmp .skip_char

.check_valid_char:
    ; Si estamos en sección de plásmido, saltar
    cmp r8, 1
    je .skip_char

    ; Validar caracteres de ADN
    cmp al, 'A'
    je .valid_char
    cmp al, 'C'
    je .valid_char
    cmp al, 'G'
    je .valid_char
    cmp al, 'T'
    je .valid_char
    cmp al, 'a'
    je .to_upper
    cmp al, 'c'
    je .to_upper
    cmp al, 'g'
    je .to_upper
    cmp al, 't'
    je .to_upper

    ; Caracter inválido
    inc qword [invalid_count]
    jmp .skip_char

.to_upper:
    sub al, 32                 ; Convertir a mayúsculas

.valid_char:
    mov [rdi], al              ; Guardar caracter válido
    inc rdi
    inc rbx                    ; Incrementar contador de caracteres limpios
    jmp .next_char

.skip_char:
    ; Mostrar mensaje de plásmido solo la primera vez
    cmp r8, 1
    jne .next_char
    cmp byte [plasmid_flag], 0
    jne .next_char
    mov byte [plasmid_flag], 1
    mov rdi, msg_plasmid
    call print_string

.next_char:
    inc rsi
    dec rcx
    jmp .process_loop

.process_done:
    mov byte [rdi], NULL       ; Terminar con NULL
    mov [clean_len], rbx       ; Guardar longitud de secuencia limpia
    mov [clean_count], rbx
    ret

; Generar k-mers y contar frecuencias
generate_kmers:
    mov rcx, [clean_len]
    mov r8, [k_value]
    sub rcx, r8                ; RCX = longitud - k
    jle .generate_done         ; Salir si no hay suficientes caracteres

    mov rsi, clean_dna         ; Puntero a secuencia limpia
    xor r9, r9                 ; Contador de posición

.generate_loop:
    ; Calcular índice para el k-mer actual
    mov rdi, rsi
    mov r10, r8                ; k_value
    xor rax, rax               ; Índice acumulado
    xor r11, r11               ; Contador interno

.index_loop:
    mov bl, [rdi]
    shl rax, 2                ; Multiplicar índice por 4 (base 4)

    ; Mapear caracter a valor
    cmp bl, 'A'
    je .a_char
    cmp bl, 'C'
    je .c_char
    cmp bl, 'G'
    je .g_char
    cmp bl, 'T'
    je .t_char

.a_char:
    ; A = 00 (0)
    jmp .next_char_index

.c_char:
    ; C = 01 (1)
    or rax, 1
    jmp .next_char_index

.g_char:
    ; G = 10 (2)
    or rax, 2
    jmp .next_char_index

.t_char:
    ; T = 11 (3)
    or rax, 3

.next_char_index:
    inc rdi
    inc r11
    cmp r11, r10
    jl .index_loop

    ; Incrementar frecuencia en tabla
    mov r11d, [freq_table + rax*4]
    inc r11d
    mov [freq_table + rax*4], r11d

    ; Siguiente posición
    inc rsi
    inc r9
    cmp r9, rcx
    jl .generate_loop

.generate_done:
    mov [total_kmers], r9
    ret

; Construir arreglo de k-mers no cero
build_kmer_array:
    mov rcx, MAX_KMERS
    xor r9, r9                ; Contador de k-mers no cero
    mov r10, [k_value]        ; Longitud k

    ; Calcular máximo índice (4^k)
    mov rax, 1
    mov rcx, r10
.calc_max_index:
    shl rax, 2
    loop .calc_max_index
    mov rcx, rax              ; RCX = max_index

    xor r11, r11              ; Índice actual

.array_loop:
    cmp r11, rcx
    jge .array_done

    ; Verificar frecuencia
    mov eax, [freq_table + r11*4]
    test eax, eax
    jz .next_index

    ; Guardar frecuencia
    mov r12, r9
    imul r12, 16
    mov [kmer_array + r12 + 12], eax

    ; Convertir índice a k-mer
    lea rdi, [kmer_array + r12] ; Buffer para cadena
    mov rsi, r11               ; Índice
    mov rdx, r10               ; k
    call decode_kmer

    ; Incrementar contador
    inc r9

.next_index:
    inc r11
    jmp .array_loop

.array_done:
    mov [num_non_zero], r9
    ret

; Ordenar k-mers con Quicksort
quicksort:
    mov r12, [num_non_zero]
    test r12, r12
    jz .done

    ; Preparar parámetros para quicksort
    xor rdi, rdi             ; left = 0
    dec r12                  ; right = num_non_zero - 1
    mov rsi, r12
    mov rdx, kmer_array
    call quicksort_recursive

.done:
    ret

quicksort_recursive:
    ; rdi = left, rsi = right, rdx = array
    cmp rdi, rsi
    jge .end_recursive

    ; Particionar
    call partition
    mov r8, rax             ; pivot_index

    ; Ordenar izquierda
    push rsi
    push r8
    mov rsi, r8
    dec rsi
    call quicksort_recursive
    pop r8
    pop rsi

    ; Ordenar derecha
    push rdi
    mov rdi, r8
    inc rdi
    call quicksort_recursive
    pop rdi

.end_recursive:
    ret

partition:
    ; rdi = left, rsi = right
    mov r9, rsi             ; pivot = right
    mov r10, rdi            ; i = left - 1
    dec r10

    mov r11, rdi            ; j = left

.partition_loop:
    cmp r11, rsi
    jge .end_partition

    ; Comparar array[j] con array[pivot]
    mov rax, r11
    imul rax, 16
    lea r12, [kmer_array + rax]

    mov rax, r9
    imul rax, 16
    lea r13, [kmer_array + rax]

    mov rdi, r12
    mov rsi, r13
    call compare_kmers
    cmp eax, 0
    jge .next_j

    ; Incrementar i e intercambiar
    inc r10
    mov rdi, r10
    mov rsi, r11
    call swap_kmers

.next_j:
    inc r11
    jmp .partition_loop

.end_partition:
    ; Intercambiar array[i+1] con array[pivot]
    inc r10
    mov rdi, r10
    mov rsi, r9
    call swap_kmers

    mov rax, r10            ; Retornar índice del pivote
    ret

; Intercambiar dos elementos en el arreglo
swap_kmers:
    ; rdi = index1, rsi = index2
    imul rdi, 16
    imul rsi, 16

    ; Intercambiar cadena (12 bytes)
    mov rax, [kmer_array + rdi]
    mov rbx, [kmer_array + rsi]
    mov [kmer_array + rsi], rax
    mov [kmer_array + rdi], rbx

    mov rax, [kmer_array + rdi + 8]
    mov rbx, [kmer_array + rsi + 8]
    mov [kmer_array + rsi + 8], rax
    mov [kmer_array + rdi + 8], rbx

    ; Intercambiar frecuencia (4 bytes)
    mov eax, [kmer_array + rdi + 12]
    mov ebx, [kmer_array + rsi + 12]
    mov [kmer_array + rsi + 12], eax
    mov [kmer_array + rdi + 12], ebx
    ret

; Buscar k-mer en tabla de frecuencias
search_kmer:
    ; rdi = puntero a cadena, rsi = k
    call encode_kmer
    cmp rax, MAX_KMERS
    jae .not_found
    
    mov ebx, [freq_table + rax*4]
    test ebx, ebx
    jz .not_found
    mov eax, ebx
    ret

.not_found:
    mov rax, -1
    ret

; ===== Funciones de utilidad =====

; Imprimir cadena (rdi = puntero a cadena)
print_string:
    push rbx
    push r12
    mov rbx, rdi
    mov rdx, 0

.count_loop:
    cmp byte [rbx], NULL
    je .count_done
    inc rdx
    inc rbx
    jmp .count_loop

.count_done:
    test rdx, rdx
    jz .print_done

    mov rax, SYS_write
    mov rsi, rdi
    mov rdi, STDOUT
    syscall

.print_done:
    pop r12
    pop rbx
    ret

; Leer entrada (rdi = buffer, rsi = tamaño máximo)
read_input:
    push rbx
    mov rbx, rdi

    ; Leer de STDIN
    mov rax, SYS_read
    mov rdi, STDIN
    mov rdx, rsi
    syscall

    ; Reemplazar LF con NULL
    mov rcx, rax
    cmp rcx, 0
    jle .read_done
    dec rcx
    mov byte [rbx + rcx], NULL

.read_done:
    pop rbx
    ret

; Convertir cadena a entero (rdi = cadena)
parse_int:
    xor rax, rax
    xor rcx, rcx

.convert_loop:
    mov cl, [rdi]
    test cl, cl
    jz .done

    cmp cl, '0'
    jb .done
    cmp cl, '9'
    ja .done

    sub cl, '0'
    imul rax, 10
    add rax, rcx
    inc rdi
    jmp .convert_loop

.done:
    ret

; Comparar cadenas (rdi = str1, rsi = str2)
strcmp:
    mov al, [rdi]
    mov bl, [rsi]
    test al, al
    jz .check_end
    cmp al, bl
    jne .not_equal
    inc rdi
    inc rsi
    jmp strcmp

.check_end:
    test bl, bl
    jnz .not_equal
    xor rax, rax
    ret

.not_equal:
    mov rax, 1
    ret

; Comparar k-mers (rdi = elem1, rsi = elem2)
compare_kmers:
    push rbx
    mov rbx, rsi

.compare_loop:
    mov al, [rdi]
    mov bl, [rsi]
    test al, al
    jz .check_end
    cmp al, bl
    jne .not_equal
    inc rdi
    inc rsi
    jmp .compare_loop

.check_end:
    test bl, bl
    jnz .not_equal
    xor eax, eax
    jmp .done

.not_equal:
    cmp al, bl
    jb .less
    mov eax, 1
    jmp .done
.less:
    mov eax, -1

.done:
    pop rbx
    ret

; Codificar k-mer a índice (rdi = cadena, rsi = k)
encode_kmer:
    xor rax, rax
    xor rcx, rcx

.encode_loop:
    cmp rcx, rsi
    jge .encode_done

    shl rax, 2
    mov bl, [rdi]

    cmp bl, 'A'
    je .next_char
    cmp bl, 'C'
    je .c_char
    cmp bl, 'G'
    je .g_char
    cmp bl, 'T'
    je .t_char

.c_char:
    or al, 1
    jmp .next_char

.g_char:
    or al, 2
    jmp .next_char

.t_char:
    or al, 3

.next_char:
    inc rdi
    inc rcx
    jmp .encode_loop

.encode_done:
    ret

; Decodificar índice a k-mer (rdi = buffer, rsi = índice, rdx = k)
decode_kmer:
    mov r8, rdx
    lea r9, [rdi + rdx]  ; Fin del buffer
    mov byte [r9], NULL  ; Terminar con NULL

.decode_loop:
    dec r9
    dec r8
    mov al, sil
    and al, 3

    cmp al, 0
    je .a_char
    cmp al, 1
    je .c_char
    cmp al, 2
    je .g_char
    cmp al, 3
    je .t_char

.a_char:
    mov byte [r9], 'A'
    jmp .next

.c_char:
    mov byte [r9], 'C'
    jmp .next

.g_char:
    mov byte [r9], 'G'
    jmp .next

.t_char:
    mov byte [r9], 'T'

.next:
    shr rsi, 2
    cmp r8, 0
    jg .decode_loop
    ret

; Imprimir entero sin signo (rdi = número)
print_uint:
    mov rax, rdi
    mov rdi, char_buffer
    mov rbx, 10
    mov rcx, 0

.convert_loop:
    xor rdx, rdx
    div rbx
    add dl, '0'
    push rdx
    inc rcx
    test rax, rax
    jnz .convert_loop

.output_loop:
    pop rax
    mov [char_buffer], al
    mov rsi, char_buffer
    mov rdx, 1
    mov rax, SYS_write
    mov rdi, STDOUT
    push rcx
    syscall
    pop rcx
    loop .output_loop
    ret
