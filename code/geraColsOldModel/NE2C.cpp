#include "NE2C.h"

NE2C::NE2C(Grafo* g, double valorSub) : valorSubtrair(valorSub), menorCustoObtido(MAIS_INFINITO){
	lixo = NULL;
	numVertices = g->getNumCustomers();
	numSatellites = g->getNumSatellites();
	vetorLabels = new vLabelNE2C[numVertices+1];
}


NE2C::~NE2C(){
	ptrLabelNE2C p, aux;
	for (int i = 1; i <= numVertices; ++i)
	{
		p = vetorLabels[i].cabeca;		
		while (p != NULL){
			aux = p;
			p = p->prox;	
			delete aux;
		}
	}
	delete [] vetorLabels;

	//apaga todos os labels que foram para o lixo (nao foram apagados pois podia haver label apontando para eles)
	p = lixo;
	while (p != NULL)
	{
		aux = p;
		p = p->prox;	
		delete aux;
	}
}


void NE2C::calculaCaminho( Grafo* g ){
	int cargaMaxima = g->getCapacVeiculosL2();
	int carga, menorIndice;
	double custo, custoAteDestino;
	ptrLabelNE2C aux, it;

	for (int i = 1; i <= numVertices; ++i)
	{
		carga = g->getCargaVertice(numSatellites + i);
		if (carga <= cargaMaxima)
		{
			aux = new LabelNE2C(i, 0, carga, g->getCustoArestaDual(0, i), NULL);
			aux->pai = NULL;
			vetorLabels[i].cabeca = vetorLabels[i].calda = vetorLabels[i].posAtual = aux;
			custoAteDestino = aux->custoDual + g->getCustoArestaDual(i, 0);
			if ((custoAteDestino < (valorSubtrair - 0.00001)) && (custoAteDestino < (menorCustoObtido - 0.00001)))
			{
				menorCustoObtido = custoAteDestino;
				menorLabelObtido = aux;
			}
		}
	}

	menorIndice = 1;
	while(menorIndice > 0)
	{
		it = vetorLabels[menorIndice].posAtual;
		vetorLabels[menorIndice].posAtual = it->prox;
		menorCustoAtual = MAIS_INFINITO;
		menorIndice = 0;

		for (int i = 1; i <= numVertices; ++i)
		{
			carga = it->cargaAcumulada + g->getCargaVertice(numSatellites + i);
			if ( (it->verticeAntecessor != i ) && ( it->ultimoVertice != i ) && ( carga <= cargaMaxima ) )
			{
				custo = it->custoDual + g->getCustoArestaDual(it->ultimoVertice, i);
				if ( !verificaLabelDominadoEDomina( i, carga, custo ) )
				{
					aux = new LabelNE2C(i, it->ultimoVertice, carga, custo, NULL);
					aux->pai = it;
					if ( vetorLabels[i].posAtual == NULL ) vetorLabels[i].posAtual = aux;
					if ( vetorLabels[i].cabeca == NULL ) vetorLabels[i].cabeca = aux;
					if ( vetorLabels[i].calda == NULL ) vetorLabels[i].calda = aux;
					else
					{
						vetorLabels[i].calda->prox = aux;
						vetorLabels[i].calda = aux;
					}
					
					
					custoAteDestino = custo + g->getCustoArestaDual(i, 0);
					if ( ( custoAteDestino < ( valorSubtrair - 0.00001 ) ) && ( custoAteDestino < ( menorCustoObtido - 0.00001 ) ) )
					{
						menorCustoObtido = custoAteDestino;
						menorLabelObtido = aux;
					}
				}
			}

			if ( ( vetorLabels[i].posAtual != NULL ) && ( vetorLabels[i].posAtual->custoDual < ( menorCustoAtual - 0.00001 ) ) )
			{
				menorIndice = i;
				menorCustoAtual = vetorLabels[i].posAtual->custoDual;
			}
		}
	}
}


bool NE2C::verificaLabelDominadoEDomina(int ult, int carga, double custo){
	ptrLabelNE2C it = vetorLabels[ult].cabeca;
	ptrLabelNE2C aux, anterior = NULL;
	bool itUpdate;
	int domCusto;

	while (it != NULL) //o novo label so sera inserido se nao houver nenhum q o domina (ou seja igual)
	{
		itUpdate = true;
		
		//Verifico se o custo deste novo label domina o label armazenado em vLabel[ult]
		if (custo < ( it->custoDual-0.00001 ) ) domCusto = 1;
		else if ( custo > (it->custoDual+0.00001 ) ) domCusto = -1;
		else domCusto = 0;

		if ( ( carga >= it->cargaAcumulada ) && ( domCusto <= 0 ) )
		{
			return true;
		}
		else if ( ( carga <= it->cargaAcumulada ) && ( domCusto >= 0 ) )
		{
			if ( vetorLabels[ult].posAtual == it ) vetorLabels[ult].posAtual = it->prox;
			if ( vetorLabels[ult].cabeca == it ) vetorLabels[ult].cabeca = it->prox;
			if ( vetorLabels[ult].calda == it ) vetorLabels[ult].calda = anterior;
			if ( anterior != NULL ) anterior->prox = it->prox;
			aux = it;
			it = it->prox;
			itUpdate = false;
			//delete aux;

			aux->prox = lixo;
			lixo = aux;
		}

		if ( itUpdate )
		{
			anterior = it;
			it = it->prox;
		}
	}
	return false;
}


vector<Rota*> NE2C::getRotaCustoMinimo(Grafo* g, int satellite, double percentual){
	Rota* r;
	vector<Rota*> rotas;
	ptrLabelNE2C aux, tmp;
	double limiteAceitacao, custoR;
	if ( menorCustoObtido < 0 ) limiteAceitacao = percentual * menorCustoObtido;
	else limiteAceitacao = ( 2 - percentual ) * menorCustoObtido;
	if ( limiteAceitacao > (valorSubtrair - 0.00001 ) ) limiteAceitacao = valorSubtrair - 0.00001;

	for (int i = 1; i <= numVertices; ++i)
	{
		aux = vetorLabels[i].cabeca;
		while (aux != NULL)
		{
			custoR = aux->custoDual + g->getCustoArestaDual(aux->ultimoVertice, 0);
			if ( custoR < limiteAceitacao )
			{
				r = new Rota();
				r->inserirVerticeFim(satellite);
				tmp = aux;
				do
				{
					r->inserirVerticeInicio(tmp->ultimoVertice+numSatellites);
					tmp = tmp->pai;
				}
				while ( tmp != NULL );
				r->inserirVerticeInicio(satellite);
				r->setCustoReduzido(custoR);
				r->setCusto(g);
				rotas.push_back(r);
			}
			aux = aux->prox;
		}
	}
	return rotas;
}
