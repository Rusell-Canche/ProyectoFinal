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
    ; Constantes del sistema
    LF              equ 10
    NULL            equ 0
    EXIT_SUCCESS    equ 0
    STDIN           equ 0
    STDOUT          equ 1
    STDERR          equ 2
    SYS_read        equ 0
    SYS_write       equ 1
    SYS_open        equ 2
    SYS_close       equ 3
    SYS_exit        equ 60
    O_RDONLY        equ 000000q
    
    ; Constantes del programa
    MAX_FILENAME    equ 256
    MAX_DNA_SIZE    equ 1000000    ; 1MB para secuencia de ADN
    MAX_KMERS       equ 100000     ; Máximo número de k-mers únicos
    MAX_KMER_LEN    equ 10
    
    ; Mensajes del programa
    welcome_msg     db LF, "=== ANALIZADOR DE ADN ===", LF
                    db "Instituto Tecnológico de Chetumal", LF
                    db "Ingeniería en Sistemas Computacionales", LF, LF, NULL
    
    filename_prompt db "Ingrese el nombre del archivo FASTA/FNA: ", NULL
    k_prompt        db "Ingrese el valor de k (4-10): ", NULL
    
    processing_msg  db LF, "Procesando archivo...", LF, NULL
    plasmid_msg     db "AVISO: El archivo contiene PLASMIDOS", LF, NULL
    invalid_char_msg db "AVISO: Caracter invalido encontrado: ", NULL
    
    results_header  db LF, "=== RESULTADOS ===", LF
                    db "K-mers encontrados y sus frecuencias:", LF, NULL
    
    search_prompt   db LF, "Buscar k-mer (ENTER para salir): ", NULL
    found_msg       db "Encontrado: aparece ", NULL
    not_found_msg   db "No encontrado", LF, NULL
    times_msg       db " veces", LF, NULL
    
    goodbye_msg     db LF, "Programa terminado.", LF, NULL
    
    ; Mensajes de error
    file_error_msg  db "Error: No se pudo abrir el archivo", LF, NULL
    k_error_msg     db "Error: k debe estar entre 4 y 10", LF, NULL
    memory_error_msg db "Error: Secuencia de ADN demasiado larga", LF, NULL
    
    ; Variables de trabajo
    filename        times MAX_FILENAME, 0
    dna_sequence    times MAX_DNA_SIZE, 0
    search_buffer   times MAX_KMER_LEN+1, 0
    temp_kmer       times MAX_KMER_LEN+1, 0
    input_buffer    times 256, 0
    
    k_value         dq 0
    dna_length      dq 0
    file_desc       dq 0
    plasmid_count   dq 0
    invalid_chars   dq 0
    
    ; Estructura para k-mers: [k-mer][frecuencia]
    ; Cada entrada ocupa MAX_KMER_LEN+1 + 8 bytes
    kmer_table      times MAX_KMERS * (MAX_KMER_LEN+1+8), 0
    kmer_count      dq 0

section .bss
    read_buffer     resb 4096

section .text
    global _start

_start:
    ; Mostrar mensaje de bienvenida
    mov rdi, welcome_msg
    call print_string
    
    ; Solicitar nombre de archivo
    mov rdi, filename_prompt
    call print_string
    
    mov rdi, filename
    mov rsi, MAX_FILENAME
    call read_string
    
    ; Solicitar valor de k
get_k_value:
    mov rdi, k_prompt
    call print_string
    
    call read_integer
    mov [k_value], rax
    
    ; Validar k (debe estar entre 4 y 10)
    cmp rax, 4
    jl invalid_k
    cmp rax, 10
    jg invalid_k
    jmp k_valid
    
invalid_k:
    mov rdi, k_error_msg
    call print_string
    jmp get_k_value
    
k_valid:
    ; Abrir archivo
    mov rax, SYS_open
    mov rdi, filename
    mov rsi, O_RDONLY
    syscall
    
    cmp rax, 0
    jl file_error
    mov [file_desc], rax
    
    ; Procesar archivo
    mov rdi, processing_msg
    call print_string
    
    call process_fasta_file
    
    ; Cerrar archivo
    mov rax, SYS_close
    mov rdi, [file_desc]
    syscall
    
    ; Generar k-mers
    call generate_kmers
    
    ; Ordenar k-mers por frecuencia
    call sort_kmers
    
    ; Mostrar resultados
    call display_results
    
    ; Modo búsqueda
    call search_mode
    
    ; Terminar programa
    mov rdi, goodbye_msg
    call print_string
    jmp exit_program

file_error:
    mov rdi, file_error_msg
    call print_string
    jmp exit_program

; Función para procesar archivo FASTA
process_fasta_file:
    push rbp
    mov rbp, rsp
    push rbx
    push r12
    push r13
    
    mov r12, 0          ; Posición en dna_sequence
    mov r13, 0          ; Flag para línea de encabezado
    
read_loop:
    ; Leer chunk del archivo
    mov rax, SYS_read
    mov rdi, [file_desc]
    mov rsi, read_buffer
    mov rdx, 4096
    syscall
    
    cmp rax, 0
    jle end_read        ; EOF o error
    
    mov rbx, 0          ; Índice en buffer
    
process_chunk:
    cmp rbx, rax
    jge read_loop
    
    mov cl, [read_buffer + rbx]
    
    ; Verificar si es línea de encabezado
    cmp cl, '>'
    je header_line
    
    ; Verificar si es salto de línea
    cmp cl, LF
    je next_char
    cmp cl, 13          ; CR
    je next_char
    
    ; Verificar si es carácter válido de ADN
    cmp cl, 'A'
    je valid_dna_char
    cmp cl, 'T'
    je valid_dna_char
    cmp cl, 'C'
    je valid_dna_char
    cmp cl, 'G'
    je valid_dna_char
    cmp cl, 'N'         ; N se considera válido pero se reporta
    je handle_n_char
    
    ; Carácter inválido
    inc qword [invalid_chars]
    call report_invalid_char
    jmp next_char
    
valid_dna_char:
    ; Verificar espacio disponible
    cmp r12, MAX_DNA_SIZE-1
    jge memory_full
    
    ; Agregar a secuencia
    mov [dna_sequence + r12], cl
    inc r12
    jmp next_char
    
handle_n_char:
    ; Tratamos N como carácter inválido pero lo reportamos
    inc qword [invalid_chars]
    jmp next_char
    
header_line:
    ; Incrementar contador de plásmidos si no es el primero
    cmp qword [plasmid_count], 0
    je first_header
    inc qword [plasmid_count]
    
first_header:
    inc qword [plasmid_count]
    
    ; Saltar hasta el final de línea
skip_header:
    inc rbx
    cmp rbx, rax
    jge read_loop
    mov cl, [read_buffer + rbx]
    cmp cl, LF
    jne skip_header
    
next_char:
    inc rbx
    jmp process_chunk
    
end_read:
    ; Terminar secuencia con NULL
    mov byte [dna_sequence + r12], NULL
    mov [dna_length], r12
    
    ; Reportar plásmidos si hay más de uno
    cmp qword [plasmid_count], 1
    jle no_plasmids
    mov rdi, plasmid_msg
    call print_string
    
no_plasmids:
    pop r13
    pop r12
    pop rbx
    pop rbp
    ret
    
memory_full:
    mov rdi, memory_error_msg
    call print_string
    jmp exit_program

report_invalid_char:
    push rax
    push rbx
    mov rdi, invalid_char_msg
    call print_string
    
    ; Mostrar el carácter inválido
    mov [temp_kmer], cl
    mov byte [temp_kmer + 1], NULL
    mov rdi, temp_kmer
    call print_string
    
    mov rdi, newline
    call print_string
    
    pop rbx
    pop rax
    ret

; Función para generar k-mers
generate_kmers:
    push rbp
    mov rbp, rsp
    push rbx
    push r12
    push r13
    push r14
    
    mov r12, 0                  ; Posición actual en ADN
    mov qword [kmer_count], 0   ; Contador de k-mers únicos
    
kmer_loop:
    ; Verificar si podemos extraer un k-mer completo
    mov rax, [dna_length]
    sub rax, [k_value]
    cmp r12, rax
    jg end_kmer_generation
    
    ; Extraer k-mer
    mov r13, 0                  ; Índice en k-mer temporal
    
extract_kmer:
    cmp r13, [k_value]
    jge kmer_extracted
    
    mov al, [dna_sequence + r12 + r13]
    mov [temp_kmer + r13], al
    inc r13
    jmp extract_kmer
    
kmer_extracted:
    ; Terminar k-mer con NULL
    mov byte [temp_kmer + r13], NULL
    
    ; Buscar si el k-mer ya existe en la tabla
    call find_or_add_kmer
    
    inc r12
    jmp kmer_loop
    
end_kmer_generation:
    pop r14
    pop r13
    pop r12
    pop rbx
    pop rbp
    ret

; Función para buscar o agregar k-mer
find_or_add_kmer:
    push rbp
    mov rbp, rsp
    push rbx
    push r12
    push r13
    
    mov r12, 0                  ; Índice en tabla de k-mers
    
search_table:
    cmp r12, [kmer_count]
    jge add_new_kmer
    
    ; Calcular posición del k-mer en tabla
    mov rax, r12
    mov rbx, MAX_KMER_LEN + 1 + 8
    mul rbx
    mov r13, rax                ; Offset en tabla
    
    ; Comparar k-mers
    mov rdi, temp_kmer
    mov rsi, kmer_table
    add rsi, r13
    call compare_strings
    
    cmp rax, 0
    je found_kmer
    
    inc r12
    jmp search_table
    
found_kmer:
    ; Incrementar frecuencia
    mov rax, r13
    add rax, MAX_KMER_LEN + 1   ; Posición de frecuencia
    inc qword [kmer_table + rax]
    jmp end_find_add
    
add_new_kmer:
    ; Verificar espacio disponible
    cmp qword [kmer_count], MAX_KMERS
    jge end_find_add
    
    ; Calcular posición para nuevo k-mer
    mov rax, [kmer_count]
    mov rbx, MAX_KMER_LEN + 1 + 8
    mul rbx
    mov r13, rax
    
    ; Copiar k-mer
    mov rdi, kmer_table
    add rdi, r13
    mov rsi, temp_kmer
    call copy_string
    
    ; Establecer frecuencia inicial en 1
    mov rax, r13
    add rax, MAX_KMER_LEN + 1
    mov qword [kmer_table + rax], 1
    
    inc qword [kmer_count]
    
end_find_add:
    pop r13
    pop r12
    pop rbx
    pop rbp
    ret

; Función de ordenamiento burbuja por frecuencia (descendente)
sort_kmers:
    push rbp
    mov rbp, rsp
    push rbx
    push r12
    push r13
    push r14
    push r15
    
    mov r12, 0                  ; i = 0
    
outer_loop:
    mov rax, [kmer_count]
    dec rax
    cmp r12, rax
    jge end_sort
    
    mov r13, 0                  ; j = 0
    
inner_loop:
    mov rax, [kmer_count]
    sub rax, r12
    dec rax
    cmp r13, rax
    jge next_outer
    
    ; Calcular posiciones de elementos a comparar
    mov rax, r13
    mov rbx, MAX_KMER_LEN + 1 + 8
    mul rbx
    mov r14, rax                ; Posición elemento j
    
    mov rax, r13
    inc rax
    mov rbx, MAX_KMER_LEN + 1 + 8
    mul rbx
    mov r15, rax                ; Posición elemento j+1
    
    ; Comparar frecuencias
    mov rax, r14
    add rax, MAX_KMER_LEN + 1
    mov rbx, [kmer_table + rax] ; Frecuencia j
    
    mov rax, r15
    add rax, MAX_KMER_LEN + 1
    mov rcx, [kmer_table + rax] ; Frecuencia j+1
    
    cmp rbx, rcx
    jge no_swap
    
    ; Intercambiar elementos completos
    call swap_kmers
    
no_swap:
    inc r13
    jmp inner_loop
    
next_outer:
    inc r12
    jmp outer_loop
    
end_sort:
    pop r15
    pop r14
    pop r13
    pop r12
    pop rbx
    pop rbp
    ret

swap_kmers:
    ; Intercambiar k-mers en posiciones r14 y r15
    push rdi
    push rsi
    push rcx
    
    mov rcx, MAX_KMER_LEN + 1 + 8   ; Tamaño total de entrada
    
swap_loop:
    cmp rcx, 0
    je end_swap
    
    mov al, [kmer_table + r14]
    mov bl, [kmer_table + r15]
    mov [kmer_table + r14], bl
    mov [kmer_table + r15], al
    
    inc r14
    inc r15
    dec rcx
    jmp swap_loop
    
end_swap:
    pop rcx
    pop rsi
    pop rdi
    ret

; Función para mostrar resultados
display_results:
    push rbp
    mov rbp, rsp
    push rbx
    push r12
    
    mov rdi, results_header
    call print_string
    
    mov r12, 0                  ; Índice actual
    
display_loop:
    cmp r12, [kmer_count]
    jge end_display
    
    ; Calcular posición del k-mer
    mov rax, r12
    mov rbx, MAX_KMER_LEN + 1 + 8
    mul rbx
    mov rbx, rax
    
    ; Mostrar k-mer
    mov rdi, kmer_table
    add rdi, rbx
    call print_string
    
    ; Mostrar espacio
    mov rdi, space_msg
    call print_string
    
    ; Mostrar frecuencia
    add rbx, MAX_KMER_LEN + 1
    mov rax, [kmer_table + rbx]
    call print_integer
    
    ; Nueva línea
    mov rdi, newline
    call print_string
    
    inc r12
    jmp display_loop
    
end_display:
    pop r12
    pop rbx
    pop rbp
    ret

; Modo de búsqueda interactivo
search_mode:
    push rbp
    mov rbp, rsp
    
search_loop:
    mov rdi, search_prompt
    call print_string
    
    mov rdi, search_buffer
    mov rsi, MAX_KMER_LEN
    call read_string
    
    ; Verificar si está vacío (salir)
    cmp byte [search_buffer], NULL
    je end_search_mode
    
    ; Buscar k-mer
    call search_kmer
    
    jmp search_loop
    
end_search_mode:
    pop rbp
    ret

search_kmer:
    push rbp
    mov rbp, rsp
    push rbx
    push r12
    
    mov r12, 0                  ; Índice en tabla
    
search_kmer_loop:
    cmp r12, [kmer_count]
    jge kmer_not_found
    
    ; Calcular posición
    mov rax, r12
    mov rbx, MAX_KMER_LEN + 1 + 8
    mul rbx
    mov rbx, rax
    
    ; Comparar
    mov rdi, search_buffer
    mov rsi, kmer_table
    add rsi, rbx
    call compare_strings
    
    cmp rax, 0
    je kmer_found_search
    
    inc r12
    jmp search_kmer_loop
    
kmer_found_search:
    mov rdi, found_msg
    call print_string
    
    ; Mostrar frecuencia
    add rbx, MAX_KMER_LEN + 1
    mov rax, [kmer_table + rbx]
    call print_integer
    
    mov rdi, times_msg
    call print_string
    jmp end_search_kmer
    
kmer_not_found:
    mov rdi, not_found_msg
    call print_string
    
end_search_kmer:
    pop r12
    pop rbx
    pop rbp
    ret

; Funciones auxiliares
print_string:
    push rbp
    mov rbp, rsp
    push rbx
    push rdx
    
    ; Contar caracteres
    mov rbx, rdi
    mov rdx, 0
    
count_loop:
    cmp byte [rbx], NULL
    je print_str
    inc rdx
    inc rbx
    jmp count_loop
    
print_str:
    mov rax, SYS_write
    mov rsi, rdi
    mov rdi, STDOUT
    syscall
    
    pop rdx
    pop rbx
    pop rbp
    ret

read_string:
    push rbp
    mov rbp, rsp
    push rbx
    push rcx
    
    mov rbx, rdi        ; Buffer destino
    mov rcx, rsi        ; Tamaño máximo
    
    mov rax, SYS_read
    mov rdi, STDIN
    mov rsi, rbx
    mov rdx, rcx
    syscall
    
    ; Remover salto de línea
    dec rax
    mov byte [rbx + rax], NULL
    
    pop rcx
    pop rbx
    pop rbp
    ret

read_integer:
    push rbp
    mov rbp, rsp
    push rbx
    
    mov rdi, input_buffer
    mov rsi, 10
    call read_string
    
    ; Convertir string a integer
    mov rdi, input_buffer
    call string_to_int
    
    pop rbx
    pop rbp
    ret

string_to_int:
    push rbp
    mov rbp, rsp
    push rbx
    
    mov rax, 0          ; Resultado
    mov rbx, 0          ; Índice
    
convert_loop:
    mov cl, [rdi + rbx]
    cmp cl, NULL
    je end_convert
    
    sub cl, '0'
    imul rax, 10
    add rax, rcx
    inc rbx
    jmp convert_loop
    
end_convert:
    pop rbx
    pop rbp
    ret

print_integer:
    push rbp
    mov rbp, rsp
    push rbx
    push rcx
    push rdx
    
    mov rbx, 10         ; Base
    mov rcx, 0          ; Contador de dígitos
    
    ; Manejar caso especial de 0
    cmp rax, 0
    jne convert_digits
    
    mov byte [input_buffer], '0'
    mov byte [input_buffer + 1], NULL
    mov rdi, input_buffer
    call print_string
    jmp end_print_int
    
convert_digits:
    mov rdx, 0
    div rbx
    add rdx, '0'
    push rdx
    inc rcx
    
    cmp rax, 0
    jne convert_digits
    
    ; Reconstruir número en buffer
    mov rbx, 0
    
build_string:
    pop rdx
    mov [input_buffer + rbx], dl
    inc rbx
    dec rcx
    cmp rcx, 0
    jne build_string
    
    mov byte [input_buffer + rbx], NULL
    mov rdi, input_buffer
    call print_string
    
end_print_int:
    pop rdx
    pop rcx
    pop rbx
    pop rbp
    ret

compare_strings:
    push rbp
    mov rbp, rsp
    push rbx
    
    mov rbx, 0
    
cmp_loop:
    mov al, [rdi + rbx]
    mov cl, [rsi + rbx]
    
    cmp al, cl
    jne strings_different
    
    cmp al, NULL
    je strings_equal
    
    inc rbx
    jmp cmp_loop
    
strings_equal:
    mov rax, 0
    jmp end_compare
    
strings_different:
    mov rax, 1
    
end_compare:
    pop rbx
    pop rbp
    ret

copy_string:
    push rbp
    mov rbp, rsp
    push rbx
    
    mov rbx, 0
    
copy_loop:
    mov al, [rsi + rbx]
    mov [rdi + rbx], al
    
    cmp al, NULL
    je end_copy
    
    inc rbx
    jmp copy_loop
    
end_copy:
    pop rbx
    pop rbp
    ret

exit_program:
    mov rax, SYS_exit
    mov rdi, EXIT_SUCCESS
    syscall

; Datos adicionales
section .data
    newline     db LF, NULL
    space_msg   db " ", NULL