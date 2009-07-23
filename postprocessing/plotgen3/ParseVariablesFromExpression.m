function [variables] = ParseVariablesFromExpression( expression )

variables =  regexp( expression, '[A-Za-z0-9_]*', 'match' );
