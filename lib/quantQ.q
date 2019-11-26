// entry point: load all code
// > cd quantQ/lib
// > q quantQ.q

{system"l ",string x}each {x where x like "quantQ_*.q"}key `:.;

-1"";
-1"Welcome to quantQ";
-1"   ____ \\ \\   ";
-1"  / ___| \\ \\  ";
-1" | |_| |  ) )   ";
-1"  \\__  | / /   ";
-1"     |_|/ /     ";
-1"";
-1"For available sub-namespaces, type key`.quantQ:";
key`.quantQ
-1"";
-1"For functions in a namespace, type .e.g \\f .quantQ.rf";
-1"";

