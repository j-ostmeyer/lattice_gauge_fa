#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*dispM[z_]:=DiagonalMatrix[Table[Total[Sin[z/2]^2],{a,z}]]-Table[Sin[a/2]Sin[b/2]Exp[I/2(a-b)],{a,z},{b,z}]*)


(* ::Input:: *)
(*dispM[{2Subscript[x, 1],2Subscript[x, 2]}]//MatrixForm*)


(* ::Input:: *)
(*Assuming[0<Subscript[x, 1]<Pi&&0<Subscript[x, 2]<Pi,Eigenvectors[dispM[{Subscript[x, 1],Subscript[x, 2]}]]//Orthogonalize//FullSimplify]*)


(* ::Input:: *)
(*Inverse[{{a,b,c},{Conjugate[b],-a,0},{Conjugate[c],0,-a}}](a^2+Abs[b]^2+Abs[c]^2)//FullSimplify*)


(* ::Input:: *)
(*Assuming[0<Subscript[x, 1]<Pi&&0<Subscript[x, 2]<Pi,PseudoInverse[dispM[{2Subscript[x, 1],2Subscript[x, 2]}]]//FullSimplify]*)


(* ::Input:: *)
(*Assuming[0<Subscript[x, 1]<Pi&&0<Subscript[x, 2]<Pi,ConjugateTranspose[Inverse[Eigenvectors[dispM[{Subscript[x, 1],Subscript[x, 2]}]]]] . dispM[{Subscript[x, 1],Subscript[x, 2]}] . Inverse[Eigenvectors[dispM[{Subscript[x, 1],Subscript[x, 2]}]]]//Simplify]*)
