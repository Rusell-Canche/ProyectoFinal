; Programa para analizar palabras de ADN en archivos FASTA/FNA
; Caracteristicas:
;   - Lee archivos FASTA/FNA, omite lineas de cabecera
;   - Genera palabras de longitud k (4-10)
;   - Guarda palabras en archivo "palabras.txt"

section .data
;--------------------------------
; Define Constants
LF equ 10                     ; line feed
NULL equ 0                    ; end of string
EXIT_SUCCESS equ 0            ; success code
STDIN equ 0                   ; standard input
STDOUT equ 1                  ; standard output
STDERR equ 2                  ; standard error
SYS_read equ 0                ; read
SYS_write equ 1               ; write
SYS_open equ 2                ; file open
SYS_close equ 3               ; file close
SYS_exit equ 60               ; terminate
SYS_creat equ 85              ; file open/create
O_RDONLY equ 000000q          ; read only
O_WRONLY equ 000001q          ; write only
O_CREAT equ 0x40
O_TRUNC equ 0x200
S_IRUSR equ 00400q            ; user read permission
S_IWUSR equ 00200q            ; user write permission
BUFF_SIZE equ 255             ; tamaño del buffer de lectura
MAX_ADN_SIZE equ 100000       ; tamaño máximo para datos de ADN (100KB)

; Mensajes y nombres de archivo
newLine db LF, NULL
header db "Analizador de palabras de ADN", LF, LF, NULL
fileName db "Prueba.txt", NULL
outFileName db "palabras.txt", NULL
promptK db "Ingrese longitud de palabra (4-10): ", NULL
errMsgOpen db "Error abriendo archivo.", LF, NULL
errMsgRead db "Error leyendo archivo.", LF, NULL
errMsgWrite db "Error escribiendo archivo.", LF, NULL
errMsgK db "Error: k debe ser entre 4 y 10.", LF, NULL
successMsg db "Palabras guardadas en palabras.txt", LF, NULL
crlf db 13, 10                ; CR + LF para compatibilidad Windows

;-------------------------------------------------------
section .bss
readBuffer resb BUFF_SIZE     ; buffer para lectura de archivo
adnBuffer resb MAX_ADN_SIZE   ; buffer para datos de ADN (sin cabeceras)
fileDesc dq 0                 ; descriptor de archivo de entrada
outFileDesc dq 0              ; descriptor de archivo de salida
kVal resb 1                   ; valor de k ingresado por usuario
adnLength resq 1              ; longitud de los datos de ADN

;++++++++++++++++++++++++++++++++++
section .text
global _start

_start:
    ; Mostrar encabezado
    mov rdi, header
    call printStr

    ; Abrir archivo de entrada
    call openInputFile
    cmp rax, 0
    jl _exitError              ; salir si hubo error

    ; Leer y procesar archivo FASTA
    call readFASTA
    cmp rax, 0
    jl _exitError              ; salir si hubo error

    ; Cerrar archivo de entrada
    mov rax, SYS_close
    mov rdi, [fileDesc]
    syscall

    ; Solicitar valor de k al usuario
    call getKValue
    cmp rax, 0
    jl _exitError              ; salir si k no es válido

    ; Generar palabras y guardar en archivo
    call generateWords
    cmp rax, 0
    jl _exitError              ; salir si hubo error

    ; Mensaje de éxito
    mov rdi, successMsg
    call printStr

    ; Salir normalmente
    mov rax, SYS_exit
    mov rdi, EXIT_SUCCESS
    syscall

_exitError:
    ; Salir con error
    mov rax, SYS_exit
    mov rdi, 1
    syscall

;-------------------------------------------------------
; Abrir archivo de entrada
; Retorna: eax = descriptor de archivo o código de error
openInputFile:
    mov rax, SYS_open
    mov rdi, fileName
    mov rsi, O_RDONLY
    syscall
    cmp rax, 0
    jl .error
    mov [fileDesc], rax
    ret
.error:
    mov rdi, errMsgOpen
    call printStr
    mov rax, -1
    ret

;-------------------------------------------------------
; Leer y procesar archivo FASTA
; Omite líneas que comienzan con '>'
readFASTA:
    xor r15, r15              ; índice para adnBuffer
    xor r12, r12              ; flag de cabecera (0 = datos, 1 = cabecera)
    
.readLoop:
    ; Leer bloque del archivo
    mov rax, SYS_read
    mov rdi, [fileDesc]
    mov rsi, readBuffer
    mov rdx, BUFF_SIZE
    syscall
    
    ; Verificar fin de archivo o error
    cmp rax, 0
    jl .readError
    je .done
    
    ; Procesar bloque leído
    mov r14, rax              ; guardar tamaño leído
    xor r13, r13              ; índice en readBuffer
    
.processByte:
    cmp r13, r14
    jge .readLoop             ; terminar bloque
    
    mov al, [readBuffer + r13]
    inc r13
    
    ; Verificar si estamos en cabecera
    test r12, r12
    jnz .headerSection
    
    ; Sección de datos de ADN
    cmp al, '>'
    je .startHeader
    cmp al, LF
    je .nextByte              ; omitir saltos de línea
    cmp al, ';'
    je .nextByte              ; omitir comentarios
    
    ; Guardar caracter válido (A,C,G,T,N)
    cmp al, 'A'
    jb .nextByte
    cmp al, 'z'
    ja .nextByte
    
    ; Convertir a mayúsculas si es necesario
    cmp al, 'a'
    jb .storeChar
    cmp al, 'z'
    ja .storeChar
    sub al, 32                ; convertir a mayúscula
    
.storeChar:
    ; Almacenar en buffer de ADN
    mov [adnBuffer + r15], al
    inc r15
    
    ; Verificar límite máximo
    cmp r15, MAX_ADN_SIZE
    jae .done                 ; detener si se alcanza el límite
    
    jmp .nextByte

.startHeader:
    mov r12, 1                ; activar flag de cabecera
    jmp .nextByte

.headerSection:
    ; En cabecera - buscar fin de línea
    cmp al, LF
    jne .nextByte
    mov r12, 0                ; desactivar cabecera

.nextByte:
    jmp .processByte

.readError:
    mov rdi, errMsgRead
    call printStr
    mov rax, -1
    ret

.done:
    ; Guardar longitud de datos de ADN
    mov [adnLength], r15
    xor rax, rax              ; retorno exitoso
    ret

;-------------------------------------------------------
; Solicitar valor de k al usuario
getKValue:
    ; Mostrar prompt
    mov rdi, promptK
    call printStr
    
    ; Leer entrada de usuario
    mov rax, SYS_read
    mov rdi, STDIN
    mov rsi, readBuffer
    mov rdx, 3                ; leer hasta 3 caracteres (k + LF)
    syscall
    
    ; Verificar si se leyó al menos 1 caracter
    cmp rax, 1
    jle .invalidK
    
    ; Convertir ASCII a entero
    mov al, [readBuffer]
    sub al, '0'
    
    ; Validar rango (4-10)
    cmp al, 4
    jb .invalidK
    cmp al, 10
    ja .invalidK
    
    ; Guardar valor de k
    mov [kVal], al
    xor rax, rax              ; retorno exitoso
    ret

.invalidK:
    mov rdi, errMsgK
    call printStr
    mov rax, -1
    ret

;-------------------------------------------------------
; Generar palabras y guardar en archivo
generateWords:
    ; Crear archivo de salida
    mov rax, SYS_creat
    mov rdi, outFileName
    mov rsi, S_IRUSR | S_IWUSR
    syscall
    cmp rax, 0
    jl .writeError
    mov [outFileDesc], rax
    
    ; Calcular parámetros
    movzx rbx, byte [kVal]    ; RBX = k (preservar)
    mov r14, [adnLength]      ; longitud ADN
    mov rax, r14
    sub rax, rbx              ; RAX = n - k
    jl .closeFile             ; si n < k, no hay palabras
    inc rax                   ; RAX = número de palabras (n - k + 1)
    mov r15, rax              ; R15 = contador de palabras
    xor r14, r14              ; índice actual en buffer ADN (R14)

    ; Generar cada palabra
.wordLoop:
    test r15, r15             ; verificar si quedan palabras
    jz .closeFile
    
    ; Escribir palabra en archivo
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    lea rsi, [adnBuffer + r14]
    mov rdx, rbx              ; longitud k
    syscall
    cmp rax, 0
    jl .writeError
    
    ; Escribir salto de línea (CR+LF para compatibilidad)
    mov rax, SYS_write
    mov rdi, [outFileDesc]
    mov rsi, crlf
    mov rdx, 2
    syscall
    cmp rax, 0
    jl .writeError
    
    ; Actualizar índices
    inc r14                   ; siguiente posición inicial
    dec r15                   ; decrementar contador de palabras
    jmp .wordLoop

.closeFile:
    ; Cerrar archivo de salida
    mov rax, SYS_close
    mov rdi, [outFileDesc]
    syscall
    xor rax, rax
    ret

.writeError:
    mov rdi, errMsgWrite
    call printStr
    mov rax, -1
    ret

;-------------------------------------------------------
; Función para imprimir cadena
; rdi = dirección de la cadena terminada en NULL
printStr:
    push rbx
    push r12
    ; Contar caracteres
    mov rbx, rdi
    mov rdx, 0
.strCountLoop:
    cmp byte [rbx], NULL
    je .strDone
    inc rdx
    inc rbx
    jmp .strCountLoop

.strDone:
    cmp rdx, 0
    je .printDone
    ; Llamar al sistema
    mov rax, SYS_write
    mov rsi, rdi
    mov rdi, STDOUT
    syscall

.printDone:
    pop r12
    pop rbx
    ret