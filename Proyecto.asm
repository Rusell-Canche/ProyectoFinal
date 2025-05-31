; Programa para analizar palabras de ADN en archivos FASTA/FNA
; Caracteristicas:
;   - Lee archivos FASTA/FNA, omite lineas de plasmidos
;   - Genera palabras de longitud k (4-10) y cuenta frecuencias
;   - Ordena palabras alfabeticamente
;   - Permite busquedas interactivas

section .data
;--------------------------------
; Define Constants
LF equ 10 ; line feed
NULL equ 0 ; end of string
TRUE equ 1
FALSE equ 0
EXIT_SUCCESS equ 0 ; success code
STDIN equ 0 ; standard input
STDOUT equ 1 ; standard output
STDERR equ 2 ; standard error
SYS_read equ 0 ; read
SYS_write equ 1 ; write
SYS_open equ 2 ; file open
SYS_close equ 3 ; file close
SYS_fork equ 57 ; fork
SYS_exit equ 60 ; terminate
SYS_creat equ 85 ; file open/create
SYS_time equ 201 ; get time
O_CREAT equ 0x40
O_TRUNC equ 0x200
O_APPEND equ 0x400
O_RDONLY equ 000000q ; read only
O_WRONLY equ 000001q ; write only
O_RDWR equ 000002q ; read and write
S_IRUSR equ 00400q
S_IWUSR equ 00200q
S_IXUSR equ 00100q
;  --------------------------------
; Define vars
BUFF_SIZE equ 255
newLine db LF, NULL
header db LF, "File Read Example."
       db LF, LF, NULL
fileName db "Prueba.txt", NULL
fileDesc  dq 0
errMsgOpen db "Error opening the file.", LF, NULL
errMsgRead db "Error reading from the file.", LF, NULL
;-------------------------------------------------------
;Memory section
section .bss
readBuffer resb BUFF_SIZE
;import file
extern fprintStr
;++++++++++++++++++++++++++++++++++
; Code section

section .text
global _start

_start:

; -----
; Display header line
mov rdi, header
call printStr

; -----
; Attempt to open file - Use system service for file open
; System Service - Open
; rax = SYS_open
; rdi = address of file name string
; rsi = attributes (i.e., read only, etc.)
; Returns:
; if error -> eax < 0
; if success -> eax = file descriptor number
; The file descriptor points to the File Control
; Block (FCB). The FCB is maintained by the OS.
; The file descriptor is used for all subsequent file
; operations (read, write, close).
openInputFile:
mov rax, SYS_open ; file open
mov rdi, fileName ; file name string
mov rsi, O_RDONLY ; read only access
syscall ; call the kernel
cmp rax, 0 ; check for success
jl errorOnOpen
mov qword [fileDesc], rax ; save descriptor

; -----
; Read from file.
; For this example, we know that the file has only 1 line.
; System Service - Read
; rax = SYS_read
; rdi = file descriptor
; rsi = address of where to place data
; rdx = count of characters to read
; Returns:
; if error -> rax < 0
; if success -> rax = count of characters actually read
mov rax, SYS_read
mov rdi, qword [fileDesc]
mov rsi, readBuffer
mov rdx, BUFF_SIZE
syscall
cmp rax, 0
jl errorOnRead

; -----
; Print the buffer.
; add the NULL for the print string
mov rsi, readBuffer
mov byte [rsi+rax], NULL
mov rdi, readBuffer
call printStr

mov rdi, newLine
call printStr


; -----
; Close the file.
; System Service - close
; rax = SYS_close
; rdi = file descriptor
mov rax, SYS_close
mov rdi, qword [fileDesc]
syscall
jmp last

; -----
; Error on open.
; note, eax contains an error code which is not used
; for this example.
errorOnOpen:
mov rdi, errMsgOpen
call printStr
jmp last
; -----
; Error on read.
; note, eax contains an error code which is not used
; for this example.
errorOnRead:
mov rdi, errMsgRead
call printStr
jmp last
last:
    mov rax, SYS_exit     ; Call code for exit
    mov rdi, EXIT_SUCCESS ; Exit program with success
    syscall


;******************************************************************
;Generic function to display a str to the screen
;str must be NULL terminated
;Method:
;   Count characters in str (excludind NULL)
;   Use syscall to output characters
global printStr
printStr:
    push rbx
    ;Count str characters
    mov rbx,rdi
    mov rdx,0
    strCountL:
        cmp byte [rbx],NULL
        je strDone
        inc rdx
        inc rbx
        jmp strCountL

    strDone:
        cmp rdx,0
        je ptrDone
        ;Call OS output str
        mov rax, SYS_write ;system cado for write()
        mov rsi, rdi ; address of chars to write
        mov rdi,STDOUT ;standard out
        syscall
    
    ptrDone:
    pop rbx
    ret