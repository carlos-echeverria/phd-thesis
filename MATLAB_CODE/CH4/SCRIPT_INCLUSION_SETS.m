% Script for performing the numerical experiments concerning the eigenvalue
% inclusion sets presented in Chapter 4 of the PH.D. thesis:
%
%  [ Echeverria - Iterative solution of discretized convection-diffusion 
%                 problems - Technische Universitaet Berlin - 2020 ]
%
%  In order to perform the experiments, it is needed to hardcode the
%  input into this file. It consist of one parameter which can be edited
%  in the section "Choose matrix to be analyzed" below.
%
% input:
%
%     Achoice: matrix to be analyzed. It accepts values 1 to 6, where:
%
%                    1 is the matrix of Example 4.15a,
%                    2 is the matrix of Example 4.15b,
%                    3 is the matrix of Example 4.16a,
%                    4 is the matrix of Example 4.16b,
%                    5 is the matrix of Example 4.17,
%                    6 is the matrix of Example 4.18,
%
% output:
%
%                   * figures for the chosen example.
%                   * data for the chosen example onto console.
%
% Written by Carlos Echeverria on May 4, 2018.
% Editted by C.E. on October  31, 2019.
% Editted by C.E. on February 24, 2020.


clearvars, close all

%% Choose matrix to be analyzed

Achoice = 1;    

%% Calculates and plots inclusion sets

inclusionsets(Achoice);
