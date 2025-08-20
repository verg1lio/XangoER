# Digital Twin de Dinâmica Veicular - Xangô e-Racing

## Introdução  
Este repositório contém o código-fonte do Digital Twin de dinâmica veicular desenvolvido pela equipe Xangô e-Racing. O software tem como objetivo auxiliar no projeto e construção de um carro de corrida elétrico do tipo FSAE Student. Por meio da modelagem e simulação de sistemas dinâmicos, o código permite realizar análises precisas de desempenho, apoiar decisões de engenharia e otimizar os parâmetros do veículo.


Nosso repósitiório completo em nosso [GitHub](https://github.com/verg1lio/XangoER)

## Sobre Dinâmica Veicular  
Dinâmica veicular é o estudo do comportamento de um veículo em movimento, abrangendo a interação entre seus diversos subsistemas. No contexto deste projeto, os principais subsistemas modelados incluem:  

- **Freios (Brake System):** Responsáveis por desacelerar o veículo de maneira controlada, garantindo segurança e eficiência.  
- **Transmissão (Drivetrain):** Conjunto de componentes que transmitem a potência do motor para as rodas, influenciando o desempenho e a eficiência energética.  
- **Suspensão: (Kinematics; Tire)** Sistema que conecta o chassi às rodas, otimizando o contato com o solo e garantindo estabilidade e conforto.  

## Estrutura do Código  
O código está estruturado em classes para representar cada subsistema dinâmico. Abaixo, um breve resumo das funcionalidades de cada classe:  

### 1. **Class Kinematics**  
Responsável por modelar a cinemática dos braços de suspensão e por caracterizar componentes como molas e amortecedores. Possibilita o acerto da dinâmica vertical do veículo para um bom equilíbrio entre segurança, desempenho e conforto.  

### 2. **Class BrakeSystem**  
Modela o sistema de freios do veículo, incluindo a força de frenagem, distribuição entre os eixos e análise de eficiência. Permite simulações para verificar o desempenho em situações críticas, como frenagens de emergência.  

### 3. **Class Drivetrain**  
Simula o funcionamento do sistema de transmissão, incluindo motor, diferencial e rodas. Aborda a conversão de potência e torque para as rodas, otimizando a entrega de força em condições variadas.  

### 4. **Class Tire**  
Representa o comportamento dos pneus por meio de modelos capazes de caracterizar a geração de forças na banda de contato, além disso permite a modelagem da geometria de direção por meio de um mecanismo de quatro barras. O modelo contribui para análises de estabilidade e controle em curvas, frenagens,acelerações e é o ponto central que integra todos os modelos dinâmicos.  

## Finalidade  
Este código é de uso interno da equipe Xangô e-Racing e foi desenvolvido exclusivamente para o projeto de um carro de corrida elétrico no âmbito da competição Formula SAE. Sua aplicação possibilita uma abordagem iterativa e fundamentada no desenvolvimento do veículo, promovendo inovação e alto desempenho.  

## Autores  
Este projeto foi desenvolvido pelos membros da equipe de Dinâmica da **Xangô e-Racing**, que contribuem com conhecimentos técnicos para o avanço do projeto.  

Lista completa de autores [[**aqui**](docs/AUTHORS.md)]

Caso precise de mais informações ou tenha interesse em contribuir, entre em contato conosco.  

---  
**Equipe Xangô e-Racing**