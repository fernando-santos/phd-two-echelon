#include "GRASP.h"

GRASP::GRASP(Grafo* G, float parametroAlfa) {
	g = G;
	alfa = parametroAlfa;
	nVertices = g->getNumCustomers();
	nSatellites = g->getNumSatellites();
	verticesForaRota = new vector<int>();
}


GRASP::~GRASP() {
	for (int i = 0; i < rotasNegativas.size(); ++i) delete rotasNegativas[i];
	delete verticesForaRota;
}


RotaG* GRASP::solucaoInicial() {
	vector<int> lrc;
	RotaG* r = new RotaG();
	int* verticeVisitado = new int[nVertices+1];
	float* estimativaCusto = new float[nVertices+1];
	float custoVisita, limite;
	float menorCusto, maiorCusto;
	int tmp, antec, ultimoVertice, tentativas = 0, capacVeiculo = g->getCapacVeiculosL2();

	//PRIMEIRA ETAPA: VERTICE QUE SAI DO DEPOSITO 0
	r->vertices.push_back(0);
	memset(verticeVisitado, 0, (nVertices+1)*sizeof(int));
	menorCusto = MAIS_INFINITO;
	maiorCusto = -MAIS_INFINITO;
	for (int i = 1; i <= nVertices; ++i)
	{
		estimativaCusto[i] = MAIS_INFINITO;
		for (int j = 0; j <= nVertices; ++j)
		{
			custoVisita = g->getCustoArestaDual( 0 , i ) + g->getCustoArestaDual( i , j );
			if (custoVisita < ( estimativaCusto[i] - 0.001 ) ) estimativaCusto[i] = custoVisita;
		}

		if (estimativaCusto[i] < ( menorCusto - 0.001 ) ) menorCusto = estimativaCusto[i];
		if (estimativaCusto[i] > ( maiorCusto + 0.001 ) ) maiorCusto = estimativaCusto[i];
	}

	//lista restrita de candidatos: elementos com custo entre [menorCusto, menorCusto + alfa*(maiorCusto-menorCusto)]
	limite = menorCusto + alfa*(maiorCusto-menorCusto);
	for (int i = 1; i <= nVertices; ++i){
		if (estimativaCusto[i] <= limite ){ //solucao deve entrar na lista restrita de candidatos
			lrc.push_back(i);
		}
	}

	//insere um vertice da lista restrica de candidatos na rota e parte para a proxima etapa da construcao	
	ultimoVertice = lrc[rand() % lrc.size()];
	verticeVisitado[ultimoVertice] = 1;
	r->vertices.push_back(ultimoVertice);
	r->custo = g->getCustoArestaDual(0, ultimoVertice);
	r->carga = g->getCargaVertice(nSatellites + ultimoVertice);
	lrc.clear();

	//SEGUNDA ETAPA: VERTICES QUE SAEM DE UM FORNECEDOR PARA OUTRO FORNECEDOR OU PARA O DEPOSITO
	do{
		//passa pelos FORNECEDORES e marca com (-1 : unreachable) aqueles que poderiam ser visitados, mas excedem a capacidade do veiculo
		for ( int i = 1; i <= nVertices; ++i )
		{
			if ( ( verticeVisitado[i] == 0 ) && ( ( r->carga + g->getCargaVertice( nSatellites + i ) ) > capacVeiculo ) ) verticeVisitado[i] = -1;
		}

		menorCusto = MAIS_INFINITO;
		maiorCusto = -MAIS_INFINITO;
		for (int i = 0; i <= nVertices; ++i)
		{
			estimativaCusto[i] = MAIS_INFINITO;
			if ( verticeVisitado[i] != 0 ) continue;

			if ( i > 0 )
			{
				for (int j = 0; j <= nVertices; ++j)
				{
					if ( verticeVisitado[j] == 0 )
					{
						custoVisita = g->getCustoArestaDual( ultimoVertice , i ) + g->getCustoArestaDual( i , j );
						if ( custoVisita < ( estimativaCusto[i] - 0.001 ) ) estimativaCusto[i] = custoVisita;
					}
				}
			}
			else
			{
				estimativaCusto[i] = g->getCustoArestaDual( ultimoVertice , 0 );
			}

			//menor e maior custos diferentes para Fornecedores e Consumidores (necessario para construir a LRC com vertices de ambos)
			if (estimativaCusto[i] < MAIS_INFINITO)
			{
				if (estimativaCusto[i] < ( menorCusto - 0.001 ) ) menorCusto = estimativaCusto[i];
				if (estimativaCusto[i] > ( maiorCusto + 0.001 ) ) maiorCusto = estimativaCusto[i];
			}
		}

		//a LRC sera composta por vertices fornecedores e consumidores, seguindo os limites de aceitacao independentes
		limite = menorCusto + alfa * ( maiorCusto - menorCusto );
		for (int i = 0; i <= nVertices; ++i)
		{
			if ( estimativaCusto[i] <= limite ) lrc.push_back(i);
		}

		//insere aleatoriamente um vertice da LRC (que pode ser consumidor ou fornecedor)
		antec = ultimoVertice;
		ultimoVertice = lrc[rand() % lrc.size()];
		verticeVisitado[ultimoVertice] = 1;
		r->vertices.push_back( ultimoVertice );
		r->custo += g->getCustoArestaDual( antec , ultimoVertice );
		r->carga += g->getCargaVertice( nSatellites + ultimoVertice );
		lrc.clear();

	}while(ultimoVertice != 0);

	//preeche o vector de vertices fora da rota que sera usado na BL insercao (por eficiencia)
	preencheVerticesForaRota( r );

	delete [] verticeVisitado;
	delete [] estimativaCusto;
	return r;	
}


void GRASP::preencheVerticesForaRota(RotaG* r){
	verticesForaRota->clear();
	int tmp, tam = r->vertices.size();
	for ( int i = 1; i <= nVertices; ++i )
	{
		tmp = 1;
		while( tmp < tam )
		{
			if (r->vertices[tmp] == i) break;
			++tmp;
		}
		if (tmp == tam) verticesForaRota->push_back(i);
	}
}


int GRASP::remocao(RotaG* r, int i){ //o valor de i passado como parametro indica onde comeca a tentar remover
	float sai, entra;
	vector<int>::iterator it = r->vertices.begin() + 1;

	while( *it != 0 )
	{
		sai = g->getCustoArestaDual( *(it-1) , *it ) + g->getCustoArestaDual( *it , *(it+1) );
		entra = g->getCustoArestaDual( *(it-1) , *(it+1) );
		if ( entra < ( sai - 0.001 ) ) //significa que o custo da solucao que remove i eh menor
		{
			r->carga -= g->getCargaVertice( nSatellites + *it );
			verticesForaRota->push_back( *it );
			r->vertices.erase( it );
			r->custo -= sai;
			r->custo += entra;
			return 1;
		}
		++it;
	}

	return 0;
}


int GRASP::insercao(RotaG* r, int i){
	float sai, entra;
	vector<int>::iterator itfr, itr = r->vertices.begin();
	int tmp, capacidadeCarga = g->getCapacVeiculosL2();
	bool firstIt = true;

	while ( ( *itr != 0 ) || ( firstIt ) )
	{
		itfr = verticesForaRota->begin();

		while ( itfr != verticesForaRota->end() )
		{
			if ( ( r->carga + g->getCargaVertice( nSatellites + *itfr ) ) <= capacidadeCarga )
			{
				sai = g->getCustoArestaDual( *itr , *(itr+1) );
				entra = g->getCustoArestaDual( *itr , *itfr ) + g->getCustoArestaDual( *itfr , *(itr+1) );
				if ( entra < ( sai - 0.001 ) )
				{
					r->carga += g->getCargaVertice( nSatellites + *itfr );
					r->vertices.insert( itr+1 , *itfr );
					verticesForaRota->erase( itfr );
					r->custo -= sai;
					r->custo += entra;
					return 1;
				}
			}
			++itfr;
		}
		++itr;
		firstIt = false;
	}

	return 0;
}


int GRASP::swap(RotaG* r){
	float entra, sai;
	int tamRota = r->vertices.size()-1;

	for ( int i = 1; i < tamRota; ++i )
	{
		for ( int j = i+1; j < tamRota; ++j )
		{
			//remove as arestas (i-1)->(i) , (i)->(i+1) , (j-1)->(j) , (j)->(j+1) [exceto se i+1 == j]
			if ( j == ( i+1 ) )
			{
				sai = g->getCustoArestaDual( r->vertices[i-1] , r->vertices[i] ) + 
						g->getCustoArestaDual( r->vertices[i] , r->vertices[i+1] ) + 
						g->getCustoArestaDual( r->vertices[j] , r->vertices[j+1] );
						
				entra = g->getCustoArestaDual( r->vertices[i-1] , r->vertices[j] ) + 
						g->getCustoArestaDual( r->vertices[j] , r->vertices[i] ) + 
						g->getCustoArestaDual( r->vertices[i] , r->vertices[j+1] );
			}
			else
			{
				sai = g->getCustoArestaDual( r->vertices[i-1] , r->vertices[i] ) + 
						g->getCustoArestaDual( r->vertices[i] , r->vertices[i+1] ) + 
						g->getCustoArestaDual( r->vertices[j-1] , r->vertices[j] ) + 
						g->getCustoArestaDual( r->vertices[j] , r->vertices[j+1] );
						
				entra = g->getCustoArestaDual( r->vertices[i-1] , r->vertices[j] ) + 
						g->getCustoArestaDual( r->vertices[j] , r->vertices[i+1] ) +
						g->getCustoArestaDual( r->vertices[j-1] , r->vertices[i] ) + 
						g->getCustoArestaDual( r->vertices[i] , r->vertices[j+1] );
			}

			if ( entra < ( sai - 0.01 ) ) //troca i e j na rota, atualiza seu custo e executa nova iteracao
			{
				int tmp = r->vertices[i];
				r->vertices[i] = r->vertices[j];
				r->vertices[j] = tmp;
				r->custo -= sai;
				r->custo += entra;
				return 1;		
			}
		}
	}

	return 0;
}


vector < Movimento > GRASP::movimentosPR(RotaG* solution, RotaG* guiding){ 	//todos os movimentos devem ser feitos na viabilidade
	vector <int> sol, guide;
	for ( int i = 1; solution->vertices[i] != 0; ++i ) sol.push_back(solution->vertices[i]);
	for ( int j = 1; guiding->vertices[j] != 0; ++j ) guide.push_back(guiding->vertices[j]);

	//obtem os possiveis movimentos para converter a solucao em guiding
	vector < Movimento > movimentos = match(sol, guide, solution->carga, 1, 1);
	return movimentos;
}


vector < Movimento > GRASP::match(vector <int>& s, vector <int>& e, int cargaS, int ajusteS, int ajusteE){
	vector < Movimento > movimentos;
	int aux, i, j, tamS = s.size(), tamE = e.size(), capacidade = g->getCapacVeiculosL2();

	//primeiro verifica todos os vertices que nao estao em solution e devem ser incluidos
	//ou aqueles que estao em solution mas nao na posicao correta (deve ser feito swap)
	for ( i = 0; i < tamE; ++i )
	{
		for ( j = 0; j < tamS; ++j )
		{
			if ( e[i] == s[j] )
			{
				//ajusta o indice para que o swap nao seja feito entre consumidores e fornecedores e vice-versa
				if ( ( ( i+ajusteE ) >= ajusteS ) && ( ( i+ajusteE ) < ( tamS+ajusteS ) ) ) 
				{
					aux = i+ajusteE;
				} else if ((i+ajusteE) >= ajusteS) {
					aux = tamS+ajusteS-1;
				} else if ( ( i+ajusteE ) <= ( tamS+ajusteS ) ) {
					aux = ajusteS;
				}

				if ( ( j != i ) && ( ( j+ajusteS ) != aux ) ) //realiza SWAP
				{
					Movimento mv ('S', j+ajusteS, aux);
					movimentos.push_back(mv);
				}
				break;
			}
		}
		if (j == tamS) //realiza a insercao do vertice que nao esta na rota
		{
			if ( cargaS + g->getCargaVertice( nSatellites + e[i] ) <= capacidade )
			{
				//ajusta o indice para que a INSERCAO nao inclua um fornecedor em um grupo de consumidores e vice-versa
				if (((i+ajusteE) >= ajusteS) && ((i+ajusteE) <= (tamS+ajusteS))) 
				{
					aux = i+ajusteE;
				} else if ((i+ajusteE) >= ajusteS) {
					aux = tamS+ajusteS;
				} else if ((i+ajusteE) <= (tamS+ajusteS)) {
					aux = ajusteS;
				}

				Movimento mv ('I', e[i], aux);
				movimentos.push_back(mv);
			}
		}
	}

	//verifica quais vertices devem ser REMOVIDOS
	for (i = 0; i < tamS; ++i)
	{
		for (j = 0; j < tamE; ++j)
		{
			if (s[i] == e[j]) break;
		}

		if (j == tamE) //significa que s[i] nao esta em guiding e deve ser removido
		{
			Movimento mv ('R', i+ajusteS);
			movimentos.push_back(mv);
		}
	}

	return movimentos;
}


RotaG* GRASP::pathRelinking(RotaG* rt1, RotaG* rt2){
	int i, j , aux, tmp;
	vector < Movimento > m;
	RotaG* melhorRotaDoCaminho = new RotaG();
	float custoMv, sai, entra, menor = MAIS_INFINITO;

	do{
		//atualiza a melhorRotaDoCaminho com as informacoes da rota cujo custo eh o menor encontrado
		if ( rt1->custo < ( menor - 0.001 ) )
		{
			melhorRotaDoCaminho->vertices.clear();
			melhorRotaDoCaminho->vertices.assign(rt1->vertices.begin(), rt1->vertices.end());
			melhorRotaDoCaminho->carga = rt1->carga;
			melhorRotaDoCaminho->custo = rt1->custo;
			menor = rt1->custo;
		}		

		m.clear();		
		m = movimentosPR(rt1, rt2);

		if ( m.size() > 0 ) 
		{			
			aux = defineIndiceMovimento(m, rt1, custoMv);
			if (m[aux].tipo == 'S') //SWAP
			{
				//realiza o movimento de swap e atualiza o custo da rota
				tmp = rt1->vertices[m[aux].arg1];
				rt1->vertices[m[aux].arg1] = rt1->vertices[m[aux].arg2];
				rt1->vertices[m[aux].arg2] = tmp;
				rt1->custo += custoMv;

			} else if (m[aux].tipo == 'I') { //INSERCAO

				//realiza o movimento de insercao e atualiza a carga e o custo da rota
				vector<int>::iterator it = rt1->vertices.begin() + m[aux].arg2;
				rt1->carga += g->getCargaVertice(nSatellites + m[aux].arg1);
				rt1->vertices.insert(it, m[aux].arg1);
				rt1->custo += custoMv;

			} else { //REMOCAO

				//realiza o movimento de remocao e atualiza a carga e o custo da rota
				vector<int>::iterator it = rt1->vertices.begin() + m[aux].arg1;
				rt1->carga -= g->getCargaVertice(nSatellites + *it);
				rt1->vertices.erase(it);
				rt1->custo += custoMv;
			}
		}
	} while ( m.size() > 0 );

	return melhorRotaDoCaminho;
}

int GRASP::defineIndiceMovimento(vector < Movimento >& mv, RotaG* sol, float& custo){
	bool movimentoViavel;
	int i, j, tmp, iMenorMv, tam = mv.size();
	float *custoMovimentos = new float[tam];
	float sai, entra, menorMv = MAIS_INFINITO;

	for ( int k = 0; k < tam; ++k )
	{
		if ( mv[k].tipo == 'S' ) //SWAP
		{
			//pega os valores dos parametros do movimento
			if (mv[k].arg1 < mv[k].arg2)
			{
				i = mv[k].arg1;
				j = mv[k].arg2;
			} else {
				i = mv[k].arg2;
				j = mv[k].arg1;
			}

			//remove as arestas (i-1)->(i) , (i)->(i+1) , (j-1)->(j) , (j)->(j+1) [exceto se i+1 == j]
			if ( j == ( i+1 ) )
			{
				sai = g->getCustoArestaDual( sol->vertices[i-1] , sol->vertices[i] ) + 
						g->getCustoArestaDual( sol->vertices[i] , sol->vertices[i+1] ) + 
						g->getCustoArestaDual( sol->vertices[j] , sol->vertices[j+1] );
						
				entra = g->getCustoArestaDual( sol->vertices[i-1] , sol->vertices[j] ) + 
						g->getCustoArestaDual( sol->vertices[j] , sol->vertices[i] ) + 
						g->getCustoArestaDual( sol->vertices[i] , sol->vertices[j+1] );
			}
			else
			{
				sai = g->getCustoArestaDual( sol->vertices[i-1] , sol->vertices[i] ) + 
						g->getCustoArestaDual( sol->vertices[i] , sol->vertices[i+1] ) + 
						g->getCustoArestaDual( sol->vertices[j-1] , sol->vertices[j] ) + 
						g->getCustoArestaDual( sol->vertices[j] , sol->vertices[j+1] );
						
				entra = g->getCustoArestaDual( sol->vertices[i-1] , sol->vertices[j] ) + 
						g->getCustoArestaDual( sol->vertices[j] , sol->vertices[i+1] ) +
						g->getCustoArestaDual( sol->vertices[j-1] , sol->vertices[i] ) + 
						g->getCustoArestaDual( sol->vertices[i] , sol->vertices[j+1] );
			}

			//armazena o custo do movimento
			custoMovimentos[k] = (entra - sai);

		}
		else if (mv[k].tipo == 'I') //INSERCAO
		{
			vector<int>::iterator it = sol->vertices.begin() + mv[k].arg2;
			sai = g->getCustoArestaDual( *(it-1) , *it );
			entra = g->getCustoArestaDual( *(it-1) , mv[k].arg1 ) + g->getCustoArestaDual( mv[k].arg1 , *it );

			//armazena o custo do movimento
			custoMovimentos[k] = (entra - sai);
		}
		else //REMOCAO
		{
			//custo de arcos com o movimento
			vector<int>::iterator it = sol->vertices.begin() + mv[k].arg1;
			sai = g->getCustoArestaDual( *(it-1) , *it ) + g->getCustoArestaDual( *it , *(it+1) );
			entra = g->getCustoArestaDual( *(it-1) , *(it+1) );

			//atualiza o custo do movimento
			custoMovimentos[k] = (entra - sai);
		}

		if (custoMovimentos[k] < menorMv){
			menorMv = custoMovimentos[k];
			iMenorMv = k;
		}
	}

	//caso o menor movimento seja Swap, abre a possibilidade de se escolher aleatoriamente ou este movimento, ou insercao ou remocao (caso existam)
	//este procedimento evita loop infinito ao escolher o swap, uma vez que este movimento nao muda a estrutura da rota quanto a presenca dos vertices
	if ( mv[iMenorMv].tipo != 'S' )
	{
		custo = custoMovimentos[iMenorMv];
	}
	else
	{
		tmp = (rand() % 2);
		if (tmp == 0) //retorna o movimento de swap mesmo
		{
			custo = custoMovimentos[iMenorMv];
		}
		else //procura pelo movimento de insercao ou remocao de menor custo, caso nao encontre, retorna o swap mesmo
		{
			menorMv = MAIS_INFINITO; tmp = -1;
			for ( int k = 0; k < tam; ++k )
			{
				if ( ( ( mv[k].tipo == 'I' ) || ( mv[k].tipo == 'R' ) )  && ( custoMovimentos[k] < menorMv ) )
				{
					menorMv = custoMovimentos[k];
					tmp = k;
				}
			}
			if ( tmp >= 0 )
			{
				custo = custoMovimentos[tmp];
				iMenorMv =  tmp;
			}
			else
			{
				custo = custoMovimentos[iMenorMv];
			}
		}
	}

	delete [] custoMovimentos;
	return iMenorMv;
}


int GRASP::run(int it_max, float valorSub, Rota* rotaCplex){
	int aux, it = 0;
	RotaG *r, *melhorPR, *guiding = NULL;
	float menorCusto = MAIS_INFINITO;

	if (rotaCplex != NULL) //vai explorar a rota obtida pelo branch and cut do cplex
	{
		//primeiro limpa as rotas obtidas da execucao da heuristica (caso existam)
		for (int i = 0; i < rotasNegativas.size(); ++i) delete rotasNegativas[i];
		rotasNegativas.clear();

		//cria uma RotaG baseando-se no objeto Rota retornado pelo cplex
		RotaG* rCplex = new RotaG();
		rCplex->carga = 0;
		rCplex->vertices.push_back(0);
		vector < int > verticesCplex = rotaCplex->getVertices();
		for ( int i = 1; verticesCplex[i] > nSatellites; ++i )
		{
			rCplex->carga += g->getCargaVertice( verticesCplex[i] );
			rCplex->vertices.push_back( verticesCplex[i] - nSatellites );
		}
		rCplex->vertices.push_back(0);
		rCplex->custo = rotaCplex->getCustoReduzido();
		insereRotaNegativa(rotasNegativas, rCplex, 'C');
		preencheVerticesForaRota(rCplex);

		//executa a busca local na solucao do cplex
		while(true){
			if ( !remocao( rCplex ) ) {
				if ( !insercao( rCplex )) {
					if ( !swap( rCplex ) ) {
						break;
					}
				}
			}
		}

		if (rCplex->custo < (rotasNegativas[0]->custo - 0.001)) insereRotaNegativa(rotasNegativas, rCplex, 'C');
		guiding = rCplex;
		menorCusto = rCplex->custo;
	}

	while(it <= it_max)
	{
		do{
			++it;
			r = solucaoInicial();
		}while(r == NULL);

		while( true )
		{
			if ( !remocao( r ) )
			{
				if ( !insercao( r ) )
				{
					if ( !swap( r ) )
					{
						break;
					}
				}
			}
		}

		if ( r->custo < ( menorCusto - 0.001 ) )
		{
			if ( guiding != NULL ) delete guiding;
			guiding = r;
			menorCusto = r->custo;
		}

		//inclue a rota obtida pela Busca Local, caso ela tenha o custo reduzido negativo
		if (r->custo < (valorSub - 0.001)) insereRotaNegativa(rotasNegativas, r, 'C');

		melhorPR = pathRelinking(r, guiding);
		if (melhorPR->custo < ( menorCusto - 0.001 ) ) 
		{ 
			//preeche o vector de vertices fora da rota que sera usado na insercao
			preencheVerticesForaRota(melhorPR);
			while(true)
			{
				if ( !remocao( melhorPR ) )
				{
					if ( !insercao( melhorPR ) )
					{
						if ( !swap(melhorPR) )
						{
							break;
						}
					}
				}
			}
			guiding = melhorPR;
			menorCusto = melhorPR->custo;
			if (guiding->custo < (valorSub - 0.001)) insereRotaNegativa(rotasNegativas, guiding, 'C');
		}

		//inclui a rota obtida pelo Path Relinking, caso ela tenha o custo reduzido negativo e seja diferente de guiding
		if (((melhorPR->custo > (guiding->custo + 0.001)) || (melhorPR->custo < (guiding->custo - 0.001))) && (melhorPR->custo < (valorSub - 0.001)))
		{
			insereRotaNegativa(rotasNegativas, melhorPR, 'P');
		}
		else if (melhorPR != guiding)
		{
			delete melhorPR;
		}

		if (r != guiding) delete r;
	}

	if ( guiding != NULL ) delete guiding;
	return rotasNegativas.size();
}

void GRASP::insereRotaNegativa(vector < RotaG* >& v, RotaG* r, char origem){
	//se a rota r for diferente das outras rotas armazenadas em v, ela sera armazenada
	bool iguais;
	int j, i = 0, tam = v.size();

	while( i < tam )
	{
		iguais = false;
		if ( r->carga == v[i]->carga ) //eh igual na carga, mas pode nao ser igual no custo
		{
			if ( ( r->custo < ( v[i]->custo + 0.01 ) ) && ( r->custo > (v[i]->custo - 0.01 ) ) ) // eh igual no custo tambem, mas pode nao ser na rota
			{
				if ( r->vertices.size() == v[i]->vertices.size() ) //eh igual em tudo, mas precisa verificar os vertices da rota
				{
					for (j = 1; j < r->vertices.size(); ++j)
					{
						if (r->vertices[j] != v[i]->vertices[j]) break;
					}
					if (j == r->vertices.size()) iguais = true;
				}
			}
		}
		if (iguais == true) break;
		++i;
	}

	if (i == tam) //significa que passou por todas as rotas de v e nenhuma delas eh igual a r, portanto, insere r
	{
		if (origem == 'C'){ //copia a rota a ser inserida, pois a original sera usada no path-relinking
			RotaG* novaRota = new RotaG();
			novaRota->vertices.assign(r->vertices.begin(), r->vertices.end());
			novaRota->carga = r->carga;
			novaRota->custo = r->custo;
			v.push_back(novaRota);
		}
		else //a rota passada como parametro sera apontada (sem fazer uma copia de r)
		{
			v.push_back(r);
		}
	}
	else
	{
		if (origem != 'C') delete r;
	}
}

Rota* GRASP::getRotaConvertida(int index, int s){
	RotaG* rG = rotasNegativas[index];
	int tamRota = rG->vertices.size()-1;

	Rota* r = new Rota();
	r->inserirVerticeFim( s );
	for ( int i = 1; i < tamRota; ++i )
	{
		r->inserirVerticeFim(nSatellites + rG->vertices[i]);
	}
	r->inserirVerticeFim( s );
	r->setCustoReduzido( rG->custo );
	r->setCusto(g);
	return r;
}
